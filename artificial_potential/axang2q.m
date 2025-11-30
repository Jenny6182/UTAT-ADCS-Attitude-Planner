function q = axang2q(axang)
% Convert axis-angle to quaternion
arguments (Input)
    axang
end

arguments (Output)
    q
end

u = axang(1:3);
theta = axang(4);

u = u / norm(u);

w = cos(theta / 2);
xyz = u * sin(theta/2);

q = [w; xyz(:)];
end