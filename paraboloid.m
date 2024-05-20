%-------------------------------------------------------%
% FEM Big Project - Antenna Structure Optimization      % 
% Eqaution of paraboloid (3D)                           %
%-------------------------------------------------------%
function z = paraboloid(x, y)
f = 2.17;
z = (x.^2 + y.^2) / (4*f);
end