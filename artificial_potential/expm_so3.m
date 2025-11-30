function q = expm_so3(u)
% Exponential map on SO(3)
% Convert a 3-vector in Lie algebra to a rotation q
arguments (Input)
    u
end

arguments (Output)
    q
end

theta = norm(u);
if theta < 1e-12
    q = [1;0;0;0];
else
    axis = u / theta;
    q = [cos(theta/2); axis*sin(theta/2)];
end
end