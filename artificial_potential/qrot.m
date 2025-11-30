function qv = qrot(q, v)
% Rotate a vector using a quaternion
arguments (Input)
    q
    v
end

arguments (Output)
    qv
end

qi = [0, v];
qv = qmult(qmult(q, qi), qinv(q));
qv = qv(2:4);
end