function q = random_quat()
% Generate a uniformly random quaternion
arguments (Input)
    
end

arguments (Output)
    q
end

x = rand([1,4]);
u = x(1);
v = x(2);
w = x(3);
q = [sqrt(u)*cos(2*pi*w);
     sqrt(1-u)*sin(2*pi*v);
     sqrt(1-u)*cos(2*pi*v);
     sqrt(u)*sin(2*pi*w)];
end