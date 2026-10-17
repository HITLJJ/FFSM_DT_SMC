function out = certified_reference_bound(z, set, geometry)
%CERTIFIED_REFERENCE_BOUND Complete certified bound for |F_i(theta)-F_i(theta0)|.
% The returned bound covers both channels for every theta in the parameter box
% and every q; z rates/accelerations are treated as fixed exact binary64 data.
% Inertia and Coriolis contributions are reported separately.

schur=schur_error_bound(set,geometry);
cor=coriolis_error_bound(z.qdot,schur);
qrdd=sym(z.qr_ddot(:),'f'); qrd=sym(z.qr_dot(:),'f');
inertia=schur.epsilon_M_exact*sum(abs(qrdd));
coriolis=cor.epsilon_C_exact*sum(abs(qrd));
total=simplify(inertia+coriolis);
out.status='CERTIFIED';
out.channel_bound_exact=[total;total];
out.channel_bound=[exact_to_upper_double(total);exact_to_upper_double(total)];
out.inertia_component_exact=inertia;
out.coriolis_component_exact=coriolis;
out.inertia_component=exact_to_upper_double(inertia);
out.coriolis_component=exact_to_upper_double(coriolis);
out.schur=schur; out.coriolis=cor;
out.scope='All theta in box, all q; fixed qdot, qr_dot and qr_ddot supplied in z.';
out.uses_samples=false;
end
