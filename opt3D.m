%-----------------------------------------------------------%
% FEM Big Project - Antenna Structure Optimization          %
% 7 Layers                                                  %
%-----------------------------------------------------------%
clc;clear

X_0 = 0.4;
A = -1;
b = 0;

options = optimoptions('fmincon', 'StepTolerance', 1e-15, 'UseParallel',true);
[X, result] = fmincon(@obj_antenna, X_0, A, b, [], [], [], [], @nonlcon, options);
[c, ~] = nonlcon(X);

l1 = 0.4; l2 = X;

%-------------- 1st layer: Hexagon (fixed) - 6 -------------%
theta_1 = linspace(0, 5/3*pi, 6)';
R_1 = l1;
x_1 = R_1 .* cos(theta_1);
y_1 = R_1 .* sin(theta_1);
z_1 = paraboloid(x_1, y_1);
x = x_1;
y = y_1;
z = z_1;

%------------- 2nd layer: dodecagon (fixed) - 12 -----------%
x_2 = zeros(12, 1);
y_2 = zeros(12, 1);
z_2 = zeros(12, 1);

% 1, 3, 5, 7, 9, 11 - On the paraboloid
Z_1 = parabola(R_1);
f = @(phi)(Z_1 + l2*sin(phi) - parabola(R_1+l2*cos(phi)));
options = optimoptions('fsolve', 'Display', 'off');
phi = fsolve(f, pi/3, options);
R_2 = R_1 + l2*cos(phi);
Z_2 = Z_1 + l2*sin(phi);

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
    f = @(X)(norm(p1-X)-l2)^2 + (norm(p2-X)-l2)^2 ...
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


%---------- 3rd layer: Icositetragon (fixed) - 24 ----------%
x_3 = zeros(24, 1);
y_3 = zeros(24, 1);
z_3 = zeros(24, 1);

index3 = 1:2:24;
[theta_2, r_2, z_2] = cart2pol(x_2, y_2, z_2);
for i=1:12
    theta = theta_2(i);
    R = r_2(i);
    Z = z_2(i);
    f = @(phi)(z + l2*sin(phi) - parabola(R+l2*cos(phi)));
    options = optimoptions('fsolve', 'Display', 'off', 'Algorithm', 'levenberg-marquardt');
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
    f = @(X)(norm(p1-X)-l2)^2 + (norm(p2-X)-l2)^2 ...
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

%------------------- 4th layer:  - 36 -----------------%
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
    Z = z_3(index3(i));
    f = @(phi)(z + l2*sin(phi) - parabola(R+l2*cos(phi)));
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
    
    f = @(X)[norm(X(1:3)-np1) - l2;
             norm(X(1:3)-np2) - l2;
             norm(X(4:6)-np3) - l2;
             norm(X(4:6)-np4) - l2;
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

%------------------- 5th layer:  - 45 -----------------%
x_5 = zeros(45, 1);
y_5 = zeros(45, 1);
z_5 = zeros(45, 1);

% 1, 4, 7, ....34
index5 = 1:5:45;
index4 = 1:4:36;
[theta_4, r_4, z_4] = cart2pol(x_4, y_4, z_4);
for i=1:9
    theta = theta_4(index4(i));
    R = r_4(index4(i));
    Z = z_4(index4(i));
    f = @(phi)(z + l2*sin(phi) - parabola(R+l2*cos(phi)));
    options = optimoptions('fsolve', 'Display', 'off', ...
                           'Algorithm', 'levenberg-marquardt');
    phi = fsolve(f, pi/3, options);
    R_5 = R + l2*cos(phi);
    Z_5 = Z + l2*sin(phi);

    z_5(index5(i)) = Z_5;
    x_5(index5(i)) = R_5 .* cos(theta);
    y_5(index5(i)) = R_5 .* sin(theta);
end

% Rest
IEN = zeros(9, 8);
for i=1:size(IEN, 1)
    IEN(i, :) = [4*i-3, 4*i-2, 4*i-2, 4*i-1, 4*i-1,4*i, ...
                 4*i, mod(4*i+1, 36)];
    np1 = [x_4(IEN(i, 1)); y_4(IEN(i, 1)); z_4(IEN(i, 1))];
    np2 = [x_4(IEN(i, 2)); y_4(IEN(i, 2)); z_4(IEN(i, 2))];
    np3 = [x_4(IEN(i, 3)); y_4(IEN(i, 3)); z_4(IEN(i, 3))];
    np4 = [x_4(IEN(i, 4)); y_4(IEN(i, 4)); z_4(IEN(i, 4))];
    np5 = [x_4(IEN(i, 5)); y_4(IEN(i, 5)); z_4(IEN(i, 5))];
    np6 = [x_4(IEN(i, 6)); y_4(IEN(i, 6)); z_4(IEN(i, 6))];
    np7 = [x_4(IEN(i, 7)); y_4(IEN(i, 7)); z_4(IEN(i, 7))];
    np8 = [x_4(IEN(i, 8)); y_4(IEN(i, 8)); z_4(IEN(i, 8))];
    np9 = [x_5(5*i-4); y_5(5*i-4); z_5(5*i-4)];
    np10 = [x_5(mod(5*i+1, 45)); y_5(mod(5*i+1, 45)); z_5(mod(5*i+1, 45))];

    f = @(X)[norm(X(1:3)-np1) - l2;
             norm(X(1:3)-np2) - l2;
             norm(X(4:6)-np3) - l2;
             norm(X(4:6)-np4) - l2;
             norm(X(7:9)-np5) - l2;
             norm(X(7:9)-np6) - l2;
             norm(X(10:12)-np7) - l2;
             norm(X(10:12)-np8) - l2;
             norm(X(1:3)-np9) - l1;
             norm(X(4:6)-X(1:3)) - l1;
             norm(X(7:9)-X(10:12)) - l1;
             norm(X(10:12)-np10) - l1];

    theta0 = linspace(0, 2*pi, 46);

    X0 = ones(12, 1);
    options = optimoptions('fsolve', 'Display', 'off', ...
                           'Algorithm', 'levenberg-marquardt');
    X = fsolve(f, X0, options);

    x_5(5*i-3) = X(1);
    y_5(5*i-3) = X(2);
    z_5(5*i-3) = X(3);
    x_5(5*i-2) = X(4);
    y_5(5*i-2) = X(5);
    z_5(5*i-2) = X(6);
    x_5(5*i-1) = X(7);
    y_5(5*i-1) = X(8);
    z_5(5*i-1) = X(9);
    x_5(5*i) = X(10);
    y_5(5*i) = X(11);
    z_5(5*i) = X(12);  
end


x = vertcat(x, x_5);
y = vertcat(y, y_5);
z = vertcat(z, z_5);
%}

zReal = paraboloid(x, y);
res = norm(zReal - z);

save pos.mat x y z l1 l2

%%
%-----------------------------------------------------------%
%-------- Plot paraboloid - polar coordinate mesh ----------%
%-----------------------------------------------------------%
n = 50;                                     % Grid density
Rmax = 3.5;
theta = linspace(0, 2*pi, n);
r = linspace(0, Rmax, n);
[R, Theta] = meshgrid(r, theta);
x0 = R .* cos(Theta);
y0 = R .* sin(Theta);
z0 = paraboloid(x0, y0);
mesh(x0, y0, z0, 'LineWidth', 1.5);hold on
colormap('bone')
scatter3(x_1, y_1, z_1, 'LineWidth', 2, 'MarkerEdgeColor', 'g'); hold on
scatter3(x_2, y_2, z_2, 'LineWidth', 2, 'MarkerEdgeColor', 'c'); hold on
scatter3(x_3, y_3, z_3, 'LineWidth', 2, 'MarkerEdgeColor', 'b'); hold on
scatter3(x_4, y_4, z_4, 'LineWidth', 2, 'MarkerEdgeColor', 'm'); hold on
scatter3(x_5, y_5, z_5, 'LineWidth', 2, 'MarkerEdgeColor', 'r'); hold on
axis equal
