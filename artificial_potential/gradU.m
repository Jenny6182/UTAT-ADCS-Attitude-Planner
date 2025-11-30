function g = gradU(q, qf, refs, obs)
% Calculate potential gradient
arguments (Input)
    q
    qf
    refs
    obs
end

arguments (Output)
    g
end

eps = 1e-6;
g = zeros(3,1);
for i=1:3
    e = zeros(3,1); e(i) = eps;
    dq_plus = expm_so3(e);
    dq_minus = expm_so3(-e);
    qp = qnormalize(qmult(dq_plus', q')); qp = qp';
    qm = qnormalize(qmult(dq_minus', q')); qm = qm';
    U_plus = potential(qp, qf, refs, obs);
    U_minus = potential(qm, qf, refs, obs);
    g(i) = (U_plus - U_minus) / (2*eps);
end
end