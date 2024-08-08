%-----------------------------------------------------------%
% FEM Big Project - Antenna Structure Optimization          %
% 7 Layers                                                  %
%-----------------------------------------------------------%
clc;clear, close all

A = [-eye(10), zeros(10,8);
    eye(10), zeros(10,8);
    zeros(8,10),-eye(8);
    zeros(8,10),eye(8)];
b= [-0.45.*ones(10,1);
    1.1.*ones(10,1);
    8.*ones(8,1);
    8.*ones(8,1)];

options = optimoptions('ga', ...
    'UseParallel', true, ...
    'MaxGenerations', 2000, ...
    'PopulationSize', 400, ...
    'FunctionTolerance', 1e-4, ...
    'ConstraintTolerance', 1e-5, ...
    'Display', 'iter');

lb = [0.2.*ones(1, 10), -5.*ones(1,8)];
ub = [1.2.*ones(1, 10), 5.*ones(1,8)];

%调用 ga 函数进行优化
[X, ~] = ga(@angle_obj_antenna_48, 18, A, b,[],[],lb,ub,[],options);

%%
l22 = X(1); l23 = X(2); l24 = X(3); l25 = X(4); l26 = X(5);
l11 = X(6); l12 = X(7); l13 = X(8); l14 = X(9); l16 = X(10);

tol1_12 = X(11)/100;
tol2_12 = X(12)/100;

tol1_23 = X(13)/100;
tol2_23 = X(14)/100;

tol1_34 = X(15)/100;
tol2_34 = X(16)/100;

tol1_56 = X(17)/100;
tol2_56 = X(18)/100;


%-------------- 1st layer: Hexagon (fixed) - 6 -------------%
l1 = l11;
theta_1 = linspace(0, 5/3*pi, 6)';
R_1 = l1;
x_1 = R_1 .* cos(theta_1);
y_1 = R_1 .* sin(theta_1);
z_1 = paraboloid(x_1, y_1);
x = x_1;
y = y_1;
z = z_1;

%------------- 2nd layer: dodecagon (fixed) - 12 -----------%
l1 = l12;
x_2 = zeros(12, 1);
y_2 = zeros(12, 1);
z_2 = zeros(12, 1);

% 1, 3, 5, 7, 9, 11 - On the paraboloid
Z_1 = parabola(R_1)-tol1_12;
f = @(phi)(Z_1 + l22*sin(phi) - parabola(R_1+l22*cos(phi))-tol2_12);
options = optimoptions('fsolve', 'Display', 'off');
phi = fsolve(f, pi/3, options);
R_2 = R_1 + l22*cos(phi);
Z_2 = Z_1 + l22*sin(phi);

index = 1:2:12;
z_2(index) = Z_2;
x_2(index) = R_2 .* cos(theta_1);
y_2(index) = R_2 .* sin(theta_1);

% 2, 4, 6, 8, 10, 12
index = 2:2:12;
IEN = [1, 2; 2, 3; 3, 4; 4, 5; 5, 6; 6, 1];
for i=index
    p1 = [x_1(IEN(i/2, 1)); y_1(IEN(i/2, 1)); z_1(IEN(i/2, 1))];
    p2 = [x_1(IEN(i/2, 2)); y_1(IEN(i/2, 2)); z_1(IEN(i/2, 2))];
    p3 = [x_2(i-1); y_2(i-1); z_2(i-1)];
    p4 = [x_2(mod(i+1, 12)); y_2(mod(i+1, 12)); z_2(mod(i+1, 12))];
    f = @(X)(norm(p1-X)-l22)^2 + (norm(p2-X)-l22)^2 ...
            + (norm(p3-X)-l1)^2 + (norm(p4-X)-l1)^2;
    X0 = ones(3, 1) .* 2;
    options = optimoptions("fsolve", "Algorithm", "levenberg-marquardt", ...
                           "StepTolerance", 1e-15, "Display", "off");
    X = fsolve(f, X0, options);
    x_2(i) = X(1);
    y_2(i) = X(2);
    z_2(i) = X(3);
end
x = vertcat(x, x_2);
y = vertcat(y, y_2);
z = vertcat(z, z_2);

if ~isreal([x,y,z])
    rms=50;
    return;
end

%---------- 3rd layer: Icositetragon (fixed) - 24 ----------%
l1 = l13;
l2 = l23;
x_3 = zeros(24, 1);
y_3 = zeros(24, 1);
z_3 = zeros(24, 1);

index3 = 1:2:24;
[theta_2, r_2, z_2] = cart2pol(x_2, y_2, z_2);
for i=1:12
    theta = theta_2(i);
    R = r_2(i);
    Z = z_2(i)-tol1_23;
    f = @(phi)(Z + l2*sin(phi) - parabola(R+l2*cos(phi))-tol2_23);
    options = optimoptions('fsolve', 'Display', 'off', ...
                           'Algorithm', 'levenberg-marquardt');
    phi = fsolve(f, pi/3, options);
    R_3 = R + l2*cos(phi);
    Z_3 = Z + l2*sin(phi);

    z_3(index3(i)) = Z_3;
    x_3(index3(i)) = R_3 .* cos(theta);
    y_3(index3(i)) = R_3 .* sin(theta);
end

% 2, 4, 6, 8, 10, 12,..., 24
index = 2:2:24;
IEN = zeros(12, 2);
for i=1:size(IEN, 1)
    IEN(i, :) = [i, mod(i,12)+1];
end
for i=index
    p1 = [x_2(IEN(i/2, 1)); y_2(IEN(i/2, 1)); z_2(IEN(i/2, 1))];
    p2 = [x_2(IEN(i/2, 2)); y_2(IEN(i/2, 2)); z_2(IEN(i/2, 2))];
    p3 = [x_3(i-1); y_3(i-1); z_3(i-1)];
    p4 = [x_3(mod(i+1, 24)); y_3(mod(i+1, 24)); z_3(mod(i+1, 24))];
    f = @(X)(norm(p1-X)-l22)^2 + (norm(p2-X)-l22)^2 ...
            + (norm(p3-X)-l1)^2 + (norm(p4-X)-l1)^2;
    X0 = ones(3, 1) .* 2;
    options = optimoptions("fsolve", "Algorithm", "levenberg-marquardt", ...
                           "StepTolerance", 1e-15, "Display", 'off');
    X = fsolve(f, X0, options);
    x_3(i) = X(1);
    y_3(i) = X(2);
    z_3(i) = X(3);
end

x = vertcat(x, x_3);
y = vertcat(y, y_3);
z = vertcat(z, z_3);

if ~isreal([x,y,z])
    rms=50;
    return;
end

%------------------- 4th layer:  - 36 -----------------%
l1 = l14;
l22 = l24;
x_4 = zeros(36, 1);
y_4 = zeros(36, 1);
z_4 = zeros(36, 1);

% 1, 4, 7, ....34
index4 = 1:3:36;
index3 = 1:2:24;
[theta_3, r_3, z_3] = cart2pol(x_3, y_3, z_3);
for i=1:12
    theta = theta_3(index3(i));
    R = r_3(index3(i));
    Z = z_3(index3(i))-tol1_34;
    f = @(phi)(Z + l2*sin(phi) - parabola(R+l2*cos(phi))-tol2_34);
    options = optimoptions('fsolve', 'Display', 'off', ...
                           'Algorithm', 'levenberg-marquardt');
    phi = fsolve(f, pi/3, options);
    R_4 = R + l2*cos(phi);
    Z_4 = Z + l2*sin(phi);

    z_4(index4(i)) = Z_4;
    x_4(index4(i)) = R_4 .* cos(theta);
    y_4(index4(i)) = R_4 .* sin(theta);
end

% Rest
IEN = zeros(12, 4);
for i=1:size(IEN, 1)
    IEN(i, :) = [2*i-1, 2*i, 2*i, mod(2*i+1, 24)];
    np1 = [x_3(IEN(i, 1)); y_3(IEN(i, 1)); z_3(IEN(i, 1))];
    np2 = [x_3(IEN(i, 2)); y_3(IEN(i, 2)); z_3(IEN(i, 2))];
    np3 = [x_3(IEN(i, 3)); y_3(IEN(i, 3)); z_3(IEN(i, 3))];
    np4 = [x_3(IEN(i, 4)); y_3(IEN(i, 4)); z_3(IEN(i, 4))];
    np5 = [x_4(3*i-2); y_4(3*i-2); z_4(3*i-2)];
    np6 = [x_4(mod(3*i+1,36)); y_4(mod(3*i+1,36)); z_4(mod(3*i+1,36))];
    
    f = @(X)[norm(X(1:3)-np1) - l22;
             norm(X(1:3)-np2) - l22;
             norm(X(4:6)-np3) - l22;
             norm(X(4:6)-np4) - l22;
             norm(X(1:3)-np5) - l1;
             norm(X(1:3)-X(4:6)) - l1;
             norm(X(4:6)-np6) - l1];
    X0 = ones(6, 1);
    options = optimoptions('fsolve', 'Display', 'off', ...
                           'Algorithm', 'levenberg-marquardt');
    X = fsolve(f, X0, options);

    x_4(3*i-1) = X(1);
    y_4(3*i-1) = X(2);
    z_4(3*i-1) = X(3);
    x_4(3*i) = X(4);
    y_4(3*i) = X(5);
    z_4(3*i) = X(6);    
end


x = vertcat(x, x_4);
y = vertcat(y, y_4);
z = vertcat(z, z_4);

if ~isreal([x,y,z])
    rms=50;
    return;
end

%------------------- 5th layer:  - 36 -----------------%
l22 = l25;
h = sqrt(l22^2-(l1/2)^2);
l22 = h;
x_5 = zeros(36, 1);
y_5 = zeros(36, 1);
z_5 = zeros(36, 1);

IEN = zeros(36, 2);
for i=1:36
    IEN(i, :) = [i, mod(i, 36)+1];
end

for i=1:size(IEN, 1)
    np1 = [x_4(IEN(i, 1)); y_4(IEN(i, 1)); z_4(IEN(i, 1))];
    np2 = [x_4(IEN(i, 2)); y_4(IEN(i, 2)); z_4(IEN(i, 2))];
    npm = 1/2.*(np1+np2);
    [theta, R, Z] = cart2pol(npm(1), npm(2), npm(3));

    f = @(phi)(Z + l22*sin(phi) - parabola(R+l22*cos(phi)));
    options = optimoptions('fsolve', 'Display', 'off', ...
                           'Algorithm', 'levenberg-marquardt');
    phi = fsolve(f, pi/3, options);
    R_5 = R + l22*cos(phi);
    Z_5 = Z + l22*sin(phi);

    z_5(i) = Z_5;
    x_5(i) = R_5 .* cos(theta);
    y_5(i) = R_5 .* sin(theta);

end


x = vertcat(x, x_5);
y = vertcat(y, y_5);
z = vertcat(z, z_5);


if ~isreal([x,y,z])
    rms=50;
    return;
end


%------------------- 6th layer:  - 48 -----------------%
l1 = l16;
l2 = l26;
x_6 = zeros(48, 1);
y_6 = zeros(48, 1);
z_6 = zeros(48, 1);

% 1, 4, 7, ....52
index6 = 1:4:48;
index5 = 1:3:36;
[theta_5, r_5, z_5] = cart2pol(x_5, y_5, z_5);
for i=1:12
    theta = theta_5(index5(i));
    R = r_5(index5(i));
    Z = z_5(index5(i));
    f = @(phi)(Z + l2*sin(phi) - parabola(R+l2*cos(phi)));
    options = optimoptions('fsolve', 'Display', 'off', ...
                           'Algorithm', 'levenberg-marquardt');
    phi = fsolve(f, pi/3, options);
    if mod(i,2) == 1
        phi = phi + tol1_56;
    else
        phi = phi + tol2_56;
    end

    R_6 = R + l2*cos(phi);
    Z_6 = Z + l2*sin(phi);

    z_6(index6(i)) = Z_6;
    x_6(index6(i)) = R_6 .* cos(theta);
    y_6(index6(i)) = R_6 .* sin(theta);
end

% Rest
IEN = zeros(12, 6);
for i=1:size(IEN, 1)
    IEN(i, :) = [3*i-2, 3*i-1, 3*i-1, 3*i, 3*i, mod(3*i+1, 36)];
    np = cell(1,size(IEN,2)+2);
    for j=1:size(IEN,2)
        np{j} = [x_5(IEN(i, j)); y_5(IEN(i, j)); z_5(IEN(i, j))];
    end
    np{end-1} = [x_6(4*i-3); y_6(4*i-3); z_6(4*i-3)];
    np{end} = [x_6(mod(4*i+1,48)); y_6(mod(4*i+1,48)); z_6(mod(4*i+1,48))];

    
    f = @(X)[norm(X(1:3)-np{1}) - l2;
             norm(X(1:3)-np{2}) - l2;
             norm(X(4:6)-np{3}) - l2;
             norm(X(4:6)-np{4}) - l2;
             norm(X(7:9)-np{5}) - l2;
             norm(X(7:9)-np{6}) - l2;
             norm(X(1:3)-np{7}) - l1;
             norm(X(7:9)-np{8}) - l1;
             norm(X(1:3)-X(4:6)) - l1;
             norm(X(7:9)-X(4:6)) - l1];
    X0 = 5 .*ones(9, 1);
    options = optimoptions('fsolve', 'Display', 'off', ...
                           'Algorithm', 'levenberg-marquardt');
    X = fsolve(f, X0, options);

    x_6(4*i-2) = X(1);
    y_6(4*i-2) = X(2);
    z_6(4*i-2) = X(3);
    x_6(4*i-1) = X(4);
    y_6(4*i-1) = X(5);
    z_6(4*i-1) = X(6);    
    x_6(4*i) = X(7);
    y_6(4*i) = X(8);
    z_6(4*i) = X(9); 
end

x = vertcat(x, x_6);
y = vertcat(y, y_6);
z = vertcat(z, z_6);

r = sqrt(x.^2 + y.^2 + z.^2);
Rm = max(r);

%%
pos = [x,y,z];
pos = [pos;[0 0 0]]; % 单位：m
num = (1:1:(size(pos,1)))'; %符号记载
points = [6, 12, 24, 36, 36, 48];
IEN = IEN_all(num, points); %可以在迭代开始的时候只生成一次，循环使用
precious_z = @(pos_xy) (pos_xy(:,1).^2+pos_xy(:,2).^2)./(4*2.17); %计算准确的抛物面句柄
rms = loss_cal(IEN, pos, precious_z);
%%
zReal = paraboloid(x, y);
res = zReal - z;
f=figure;
stem(1:length(z), res, 'LineWidth', 1.5);
saveas(f, 'data/res', 'fig');

save data/pos.mat x y z res rms

%%
%-----------------------------------------------------------%
%-------- Plot paraboloid - polar coordinate mesh ----------%
%-----------------------------------------------------------%
f=figure;
n = 50;                                     % Grid density
Rmax = Rm;
theta = linspace(0, 2*pi, n);
r = linspace(0, Rmax, n);
[R, Theta] = meshgrid(r, theta);
x0 = R .* cos(Theta);
y0 = R .* sin(Theta);
z0 = paraboloid(x0, y0);
mesh(x0, y0, z0, 'LineWidth', 1.5);hold on
scatter3(x_1, y_1, z_1, 'LineWidth', 2, 'MarkerEdgeColor', 'c'); hold on
scatter3(x_2, y_2, z_2, 'LineWidth', 2, 'MarkerEdgeColor', 'b'); hold on
scatter3(x_3, y_3, z_3, 'LineWidth', 2, 'MarkerEdgeColor', 'g'); hold on
scatter3(x_4, y_4, z_4, 'LineWidth', 2, 'MarkerEdgeColor', 'm'); hold on
scatter3(x_5, y_5, z_5, 'LineWidth', 2, 'MarkerEdgeColor', 'r'); hold on
scatter3(x_6, y_6, z_6, 'LineWidth', 2, 'MarkerEdgeColor', 'k'); hold on
colorbar
axis equal
saveas(f, 'data/antenna', 'fig')

