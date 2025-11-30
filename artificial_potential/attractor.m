function U = attractor(k, qr, qt)
% Attractive potential function between two quaternions using logmap in
% SO(3)
arguments (Input)
    k
    qr
    qt
end

arguments (Output)
    U
end

U = k * norm(qlogmap(qr, qt))^2;
end