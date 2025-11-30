function xi = q2xi(q)
% Convert a quaternion into an axis-angle vector;
arguments (Input)
    q
end

arguments (Output)
    xi
end

w = q(1);
v = q(2:4);  % vector part

% Ensure unit quaternion
qnorm = norm(q);
if abs(qnorm-1) > 1e-12
    q = q / qnorm;
    w = q(1);
    v = q(2:4);
end

vnorm = norm(v);

% Handle small angles
if vnorm < 1e-12
    xi = [0; 0; 0];
    return;
end

% Rotation angle
theta = 2 * atan2(vnorm, w);

% Rotation axis
axis = v / vnorm;

% Axis-angle vector
xi = theta * axis;
end