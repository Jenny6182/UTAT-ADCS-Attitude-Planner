function qn = qnormalize(q)
% Normalize a quaternion to a unit quaternion
arguments (Input)
    q
end

arguments (Output)
    qn
end

n = norm(q);
if n < 1e-12
    qn = [1; 0; 0; 0];
else
    qn = q / n;
end
end