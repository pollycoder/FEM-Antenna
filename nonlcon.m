function [c, ceq] = nonlcon(X)
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
                           "StepTolerance", 1e-15, "Display", "off");
    X = fsolve(f, X0, options);
    x_3(i) = X(1);
    y_3(i) = X(2);
    z_3(i) = X(3);
end

x = vertcat(x, x_3);
y = vertcat(y, y_3);
z = vertcat(z, z_3);

r_3 = sqrt(x_3.^2 + y_3.^2 + z_3.^2);


c = [R_2 .* ones(size(r_2)) - r_2; R_3 .* 0.9 .* ones(size(r_3)) - r_3];
ceq = [];

end