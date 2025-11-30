function qi = qinv(q)
% Inverse of a quaternion
% Assumes Hamilton, scalar-first convention
arguments (Input)
    q
end

arguments (Output)
    qi
end

qi = [q(1), -q(2), -q(3), -q(4)];
end