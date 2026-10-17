function out = coriolis_error_bound(qdot, schurCertificate)
%CORIOLIS_ERROR_BOUND Certified ||C(theta)-C(theta0)||_2 bound at qdot.
% qdot is the fixed physical joint velocity [rad/s]. The bound uses the
% certified global derivative errors of Mr and exact rational arithmetic.

assert(strcmp(schurCertificate.status,'CERTIFIED'),'FFSM:CertificateRequired', ...
    'A certified Schur derivative envelope is required.');
qdot=sym(qdot(:),'f');
assert(numel(qdot)==2,'FFSM:Dimension','qdot must contain two rates.');
g=schurCertificate.mass_derivative_error_exact(:);
E=repmat(sym(0),2,2);
for i=1:2
    for j=1:2
        for k=1:2
            E(i,j)=E(i,j)+sym(1)/2*(g(k)+g(j)+g(i))*abs(qdot(k));
        end
    end
end
% ||X||_2 <= ||X||_F <= sum_ij |X_ij|; the latter avoids square roots and
% remains exact rational for exact dyadic input data.
epsC=sum(E,'all');
out.status='CERTIFIED';
out.element_bounds_exact=simplify(E);
out.epsilon_C_exact=simplify(epsC);
out.epsilon_C=exact_to_upper_double(out.epsilon_C_exact);
out.proof_basis='Christoffel entry bounds plus spectral <= Frobenius <= entrywise 1-norm.';
end
