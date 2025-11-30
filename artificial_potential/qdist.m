function d = qdist(q1, q2)
% Calculate the rotational distance between two quaternions
arguments (Input)
    q1
    q2
end

arguments (Output)
    d
end
i = sum(dot(q1, q2));
d = 2*acos(i);
end