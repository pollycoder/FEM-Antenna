%-------------------------------------------------------%
% FEM Big Project - Antenna Structure Optimization      % 
% Eqaution of parabola (1D)                             %
%-------------------------------------------------------%
function z = parabola(r)
f = 2.17;
z = r.^2 / (4*f);
end