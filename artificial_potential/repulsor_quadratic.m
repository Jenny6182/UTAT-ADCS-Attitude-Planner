function [U,Up] = repulsor_quadratic(u_r, u_o, t_min, k)
% Repulsive potential source with quadratic formulation
arguments (Input)
    u_r
    u_o
    t_min
    k
end

arguments (Output)
    U
    Up
end

d = dist_angle(u_r, u_o);
if d < t_min
    U = 0.5 * k * (t_min - d)^2;
    Up = -k * (t_min - d);
else
    U = 0.0;
    Up = 0.0;
end

end