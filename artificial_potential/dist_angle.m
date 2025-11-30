function t = dist_angle(u1, u2)
% Angular distance between two vectors
arguments (Input)
    u1
    u2
end

arguments (Output)
    t
end

dprod = dot(u1, u2) / (norm(u1) * norm(u2));
dprod = max(-1, min(1, dprod));
crossnorm = norm(cross(u1, u2));
t = atan2(crossnorm, dprod);
end