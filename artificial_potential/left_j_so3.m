function J = left_j_so3(xi)
% Left Jacobian of SO(3)
arguments (Input)
    xi
end

arguments (Output)
    J
end

theta = norm(xi);
if theta < 1e-8
    J = eye(3) - 0.5 * skew(xi) + (1/12) * (skew(xi)^2);
    return
end
A = (1 - cos(theta)) / (theta^2);
B = (theta - sin(theta)) / (theta^3);
K = skew(xi);
J = eye(3) - A * K + B * (K*K);
end