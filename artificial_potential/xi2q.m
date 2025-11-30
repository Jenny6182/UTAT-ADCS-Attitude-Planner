function q = xi2q(xi)
% Convert axis-angle vector to quaternion
arguments (Input)
    xi
end

arguments (Output)
    q
end

t = norm(xi);

if t < 1e-12
    q = [1; 0; 0; 0];
else
    axis = xi / t;
    half = t / 2;
    q = [cos(half);
         axis(1)*sin(half);
         axis(2)*sin(half);
         axis(3)*sin(half)];
end

end