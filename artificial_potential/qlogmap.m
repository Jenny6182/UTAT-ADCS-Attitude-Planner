function l = qlogmap(q1, q2)
% Log map in SO(3) from quaternions
arguments (Input)
    q1
    q2
end

arguments (Output)
    l
end
% relative rotation
q = qmult(qinv(q1), q2);

% shortest rotaiton
if q(1) < 0
    q = -q;
end

% clip numerical errors
q = max(-1, min(q, 1.0));

t = 2 * acos(q(1));

s = sqrt(1 - q(1)*q(1));
u = [q(2), q(3), q(4)] ./ s;

l = t*u;

end