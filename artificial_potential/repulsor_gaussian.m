function [U,Up] = repulsor_gaussian(u_r, u_o, k, s)
% Repuslive potential source with gaussian formulation
arguments (Input)
    u_r
    u_o
    k
    s
end

arguments (Output)
    U
    Up
end

d = dist_angle(u_r, u_o);

U = k * exp(-0.5 * (d./s).^2);
Up = -k * (d./(s^2)) .* exp(-0.5 * (d./s).^2);

end