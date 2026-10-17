function out = schur_error_bound(set, geometry)
%SCHUR_ERROR_BOUND Certified uniform Schur-complement inertia error bound.
% Returns ||Mr(q,theta)-Mr(q,theta0)||_2 <= epsilon_M for every q and every
% theta in the box. All certificate arithmetic is exact rational symbolic
% arithmetic; displayed doubles are outward-enclosing conversions.

assert(license('test','Symbolic_Toolbox'),'FFSM:SymbolicUnavailable', ...
    'Symbolic Math Toolbox is required for CERTIFIED arithmetic.');
c=analytic_block_norms(geometry,true);
t0=sym(set.nominal,'f'); lo=sym(set.lower,'f'); hi=sym(set.upper,'f');
r=sym(set.halfwidth,'f');
dotb=@(x,y) sum(x(:).*y(:));

A0=dotb(t0,c.A); B0=dotb(t0,c.B); D0=dotb(t0,c.D);
Amax=dotb(hi,c.A); Bmax=dotb(hi,c.B); Dmax=dotb(hi,c.D); %#ok<NASGU>
dA=dotb(r,c.A); dB=dotb(r,c.B); dD=dotb(r,c.D);

% Base mass and base COM inertia alone give A >= diag(m0,m0,J0).
k0=1/min(t0(1),t0(4));
kmax=1/min(lo(1),lo(4));
dK=k0*dA*kmax;
epsS=dB*kmax*Bmax + B0*dK*Bmax + B0*k0*dB;
epsM=dD+epsS;

g=repmat(sym(0),2,1);
for k=1:2
    Ak0=dotb(t0,c.Aq(:,k)); Bk0=dotb(t0,c.Bq(:,k));
    Akmax=dotb(hi,c.Aq(:,k)); Bkmax=dotb(hi,c.Bq(:,k));
    dAk=dotb(r,c.Aq(:,k)); dBk=dotb(r,c.Bq(:,k)); dDk=dotb(r,c.Dq(:,k));
    termLeft=dBk*kmax*Bmax + Bk0*dK*Bmax + Bk0*k0*dB;
    termRight=dB*kmax*Bkmax + B0*dK*Bkmax + B0*k0*dBk;
    termFive=dB*kmax*Akmax*kmax*Bmax ...
        +B0*dK*Akmax*kmax*Bmax ...
        +B0*k0*dAk*kmax*Bmax ...
        +B0*k0*Ak0*dK*Bmax ...
        +B0*k0*Ak0*k0*dB;
    g(k)=dDk+termLeft+termRight+termFive;
end

out.status='CERTIFIED';
out.proof_basis=['Exact-rational global Jacobian block bounds, base-only ', ...
    'Loewner lower bound, inverse perturbation identity, and norm submultiplicativity.'];
out.epsilon_A_exact=simplify(dA); out.epsilon_B_exact=simplify(dB);
out.epsilon_D_exact=simplify(dD); out.inverse_nominal_exact=simplify(k0);
out.inverse_global_exact=simplify(kmax); out.inverse_error_exact=simplify(dK);
out.epsilon_M_exact=simplify(epsM);
out.mass_derivative_error_exact=simplify(g);
out.epsilon_M=exact_to_upper_double(out.epsilon_M_exact);
out.mass_derivative_error=arrayfun(@exact_to_upper_double,out.mass_derivative_error_exact);
out.nominal_block_bounds=double([A0;B0;D0]);
out.global_block_bounds=double([Amax;Bmax;dotb(hi,c.D)]);
end
