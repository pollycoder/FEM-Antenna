%-------------------------------------------------------------------%
% Big Project - Flexible Antenna Structure Optimization             %
%-------------------------------------------------------------------%
clc;clear

% Bounds
l_lb = 0.4;
l_ub = 0.6;
theta_lb = 0;
theta_ub = pi/3;
lb = [l_lb, theta_lb];
ub = [l_ub, theta_ub];

% Initial values
X0 = [0.45, pi/30];

% Optimization
[X, res] = fmincon(@obj_linear_approx, X0, [],[],[],[],lb,ub);

%-------------------------------------------------------------------%
%-------------------- Calculation the Variables --------------------%
%-------------------------------------------------------------------%
l = X(1); theta = X(2);
rmax = 3;
zmax = parabola(rmax);

% Radius and vector list
r = zeros(100, 1);
c = zeros(100, 1);
z = zeros(100, 1);

i0 = 1;
theta_temp = 0;
while(true)
    theta_temp = theta_temp + theta;
    c(i0+1) = c(i0) + l * exp(theta_temp*1i);
    r(i0+1) = real(c(i0+1));
    z(i0+1) = imag(c(i0+1));
    theta_c = angle(c(i0+1));
    i0 = i0 + 1;
    if r(i0) >= rmax 
        break;
    elseif theta_c >= atan(zmax/rmax)
        break;
    elseif i0 == length(r)
        break;
    end
end


r = r(r>0);
r = [0; r];
r = r(r<3);
z = z(1:length(r));
z_real = parabola(r);

res_vec = z_real - z;
res = sum(abs(res_vec))/length(res_vec);


%-------------------------------------------------------------------%
%------------------------------- Plot ------------------------------%
%-------------------------------------------------------------------%

% Baseline
r_base = linspace(0, rmax, 100);
z_base = parabola(r_base);
plot(r_base, z_base, 'r-', 'LineWidth', 1.5);hold on

% Samples
plot(r, z_real, 'ro', 'LineWidth', 1.5); hold on

% Approximated curve
plot(r, z, 'bv-', 'LineWidth', 1.5);hold on
legend('Baseline', 'Sample points', 'Approximated points')
axis equal


%-------------------------------------------------------------------%
%-------------------------- Functions ------------------------------%
%-------------------------------------------------------------------%

% Objective function
function res = obj_linear_approx(X)
l = X(1); theta = X(2);
rmax = 3;
zmax = parabola(rmax);

% Radius and complex number list
r = zeros(100, 1);
c = zeros(100, 1);
z = zeros(100, 1);

i0 = 1;
theta_temp = 0;
while(true)
    theta_temp = theta_temp + theta;
    c(i0+1) = c(i0) + l * exp(theta_temp*1i);
    r(i0+1) = real(c(i0+1));
    z(i0+1) = imag(c(i0+1));
    theta_c = angle(c(i0+1));
    i0 = i0 + 1;
    if r(i0) >= rmax 
        break;
    elseif theta_c >= atan(zmax/rmax)
        break;
    elseif i0 == length(r)
        break;
    end
end

r = r(r>0);
r = [0; r];
r = r(r<3);
z = z(1:length(r));
z_real = parabola(r);

res_vec = z_real - z;
res = sum(abs(res_vec))/length(res_vec);
end