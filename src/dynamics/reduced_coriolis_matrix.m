function C = reduced_coriolis_matrix(q, qdot, p)
%REDUCED_CORIOLIS_MATRIX Christoffel matrix for the zero-momentum metric.
% Complex-step derivatives avoid subtractive cancellation and do not invoke
% symbolic operations inside an ODE step.

q = q(:); qdot = qdot(:);
assert(numel(q) == 2 && numel(qdot) == 2, 'FFSM:Dimension', ...
    'q and qdot must each contain two elements.');
step = 1e-20;
dM = zeros(2,2,2);
for k = 1:2
    qc = complex(q);
    qc(k) = qc(k) + 1i*step;
    dM(:,:,k) = imag(reduced_mass_matrix(qc, p)) / step;
end

C = zeros(2,2);
for i = 1:2
    for j = 1:2
        for k = 1:2
            C(i,j) = C(i,j) + 0.5 * (dM(i,j,k) + ...
                dM(i,k,j) - dM(j,k,i)) * qdot(k);
        end
    end
end
end
