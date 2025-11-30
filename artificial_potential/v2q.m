function q = v2q(v,u)
% Map a vector v to a quaternion that rotates u to v
arguments (Input)
    v
    u
end

arguments (Output)
    q
end

% Map v on sphere to a quaternion that rotates u to v
v = v(:)/norm(v);

% If v == b or v == -b handle separately
if norm(v - u) < 1e-6
    q = [1;0;0;0];
    return;
elseif norm(v + u) < 1e-6
    q = [0;0;1;0];     % 180 deg around y axis
    return;
end

axis = cross(u, v);
angle = acos(dot(u, v));
axis = axis ./ norm(axis);
q = axang2q([axis', angle])';
end