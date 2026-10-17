function out=certified_hessian_bound(z,set,geometry,precomputedCoeff)
%CERTIFIED_HESSIAN_BOUND Strict global Taylor-Hessian torque-error bound.
% For fixed rates and accelerations in z, this routine bounds every entry of
% d^2F_i/dtheta^2 over the complete positive parameter box and every q.
% Exact symbolic rational arithmetic and outward-checked binary64 conversion
% are used. No samples, optimizers, or unverified interval endpoints enter.

tStart=tic;
out=failed_result(z,set);
try
    assert(license('test','Symbolic_Toolbox'),'FFSM:SymbolicUnavailable', ...
        'Symbolic Math Toolbox is required for exact certificate arithmetic.');
    assert(all(set.lower>0),'FFSM:InvalidParameterSet','The complete box must be positive.');
    if nargin<4, c=fixed_state_block_norms(z.q,geometry); else, c=precomputedCoeff; end
    lo=sym(set.lower(:),'f'); hi=sym(set.upper(:),'f');
    theta0=sym(set.nominal(:),'f'); r=sym(set.halfwidth(:),'f');
    qdot=sym(z.qdot(:),'f'); qrd=sym(z.qr_dot(:),'f');
    qrdd=sym(z.qr_ddot(:),'f');

    % All link rotational inertias add positive e3*e3' base blocks.
    kappa=1/min(lo(1),sum(lo(4:6)));
    B=sum(hi.*c.B); Bq=(hi.'*c.Bq).'; Aq=(hi.'*c.Aq).';
    B0=sum(theta0.*c.B); Bq0=(theta0.'*c.Bq).'; Aq0=(theta0.'*c.Aq).';
    kappa0=1/min(theta0(1),sum(theta0(4:6)));
    gm=sym(zeros(6,1)); gmq=sym(zeros(6,2));
    um=sym(zeros(6,6)); umq=sym(zeros(6,6,2));
    for j=1:6
        gm(j)=c.D(j)+first_bound(j,0,c,kappa0,B0,Bq0,Aq0);
        for ell=1:2
            gmq(j,ell)=c.Dq(j,ell)+first_bound(j,ell,c,kappa0,B0,Bq0,Aq0);
        end
        for kk=1:6
            um(j,kk)=second_bound(j,kk,0,c,kappa,B,Bq,Aq);
            for ell=1:2
                umq(j,kk,ell)=second_bound(j,kk,ell,c,kappa,B,Bq,Aq);
            end
        end
    end

    Ec=sym(zeros(2,2,6)); Uc=sym(zeros(2,2,6,6));
    for a=1:2
        for b=1:2
            for j=1:6
                for ell=1:2
                    Ec(a,b,j)=Ec(a,b,j)+sym(1)/2*(gmq(j,ell)+ ...
                        gmq(j,b)+gmq(j,a))*abs(qdot(ell));
                    for kk=1:6
                        Uc(a,b,j,kk)=Uc(a,b,j,kk)+sym(1)/2*( ...
                            umq(j,kk,ell)+umq(j,kk,b)+umq(j,kk,a))*abs(qdot(ell));
                    end
                end
            end
        end
    end

    G=sym(zeros(2,6)); U=sym(zeros(6,6,2));
    for channel=1:2
        for j=1:6
            G(channel,j)=gm(j)*sum(abs(qrdd));
            for col=1:2
                G(channel,j)=G(channel,j)+Ec(channel,col,j)*abs(qrd(col));
            end
            for kk=1:6
                U(j,kk,channel)=um(j,kk)*sum(abs(qrdd));
                for col=1:2
                    U(j,kk,channel)=U(j,kk,channel)+ ...
                        Uc(channel,col,j,kk)*abs(qrd(col));
                end
            end
        end
    end

    linear=sym(zeros(2,1)); element=sym(zeros(2,1));
    rawSpectral=sym(zeros(2,1)); normalizedSpectral=sym(zeros(2,1));
    Lraw=sym(zeros(2,1)); Lnormalized=sym(zeros(2,1));
    rxi=r./theta0;
    for channel=1:2
        linear(channel)=sum(G(channel,:).'.*r);
        element(channel)=sym(1)/2*sum(sum(U(:,:,channel).*(r*r.')));
        Lraw(channel)=sqrt(sum(U(:,:,channel).^2,'all'));
        rawSpectral(channel)=sym(1)/2*Lraw(channel)*sum(r.^2);
        Uxi=U(:,:,channel).*(theta0*theta0.');
        Lnormalized(channel)=sqrt(sum(Uxi.^2,'all'));
        normalizedSpectral(channel)=sym(1)/2*Lnormalized(channel)*sum(rxi.^2);
    end
    remainder=sym(zeros(2,1)); chosen=cell(2,1);
    for channel=1:2
        values=[element(channel);rawSpectral(channel);normalizedSpectral(channel)];
        [remainder(channel),index]=exact_min(values);
        names={'ELEMENTWISE','RAW_SPECTRAL','NORMALIZED_SPECTRAL'};
        chosen{channel}=names{index};
    end

    out.status='CERTIFIED';
    out.bound=upper_double(linear+remainder);
    out.linear_part=upper_double(linear);
    out.remainder_part=upper_double(remainder);
    out.hessian_upper_bound=upper_double(U);
    out.spectral_hessian_upper_bound=upper_double(Lraw);
    out.normalized_spectral_hessian_upper_bound=upper_double(Lnormalized);
    out.elementwise_remainder=upper_double(element);
    out.raw_spectral_remainder=upper_double(rawSpectral);
    out.normalized_spectral_remainder=upper_double(normalizedSpectral);
    out.chosen_remainder_method=chosen;
    out.parameter_gradient_upper_bound=upper_double(G);
    out.parameter_scale=double(set.nominal(:));
    out.normalized_radius=double(set.halfwidth(:)./set.nominal(:));
    out.bound_exact=simplify(linear+remainder);
    out.linear_part_exact=simplify(linear);
    out.remainder_part_exact=simplify(remainder);
    out.hessian_upper_bound_exact=simplify(U);
    out.parameter_box=set;
    out.state=z;
    out.proof_method=['Exact-rational fixed-state trigonometric interval enclosure; ', ...
        'Schur product-word differentiation; base Loewner inverse bound; Christoffel propagation; ', ...
        'elementwise and Frobenius Taylor remainder bounds; checked outward conversion.'];
    out.failure_reason='';
    out.uses_samples=false;
    out.subbox_count=1;
catch exception
    out.status='FAILED';
    out.failure_reason=sprintf('%s: %s',exception.identifier,exception.message);
end
out.runtime=toc(tStart);
end

function value=first_bound(j,ell,c,kappa,B,Bq,Aq)
fB=factor(B,select(Bq,ell)); fBj=factor(c.B(j),select(c.Bq(j,:),ell));
fAj=factor(c.A(j),select(c.Aq(j,:),ell)); fK=factor(kappa,select(kappa^2*Aq,ell));
value=word([fBj,fK,fB])+word([fB,fK,fBj])+word([fB,fK,fAj,fK,fB]);
if ell>0, value=wordq([fBj,fK,fB])+wordq([fB,fK,fBj])+ ...
        wordq([fB,fK,fAj,fK,fB]); end
end

function value=second_bound(j,k,ell,c,kappa,B,Bq,Aq)
fB=factor(B,select(Bq,ell)); fJ=factor(c.B(j),select(c.Bq(j,:),ell));
fKp=factor(c.B(k),select(c.Bq(k,:),ell));
fAj=factor(c.A(j),select(c.Aq(j,:),ell));
fAk=factor(c.A(k),select(c.Aq(k,:),ell));
fInv=factor(kappa,select(kappa^2*Aq,ell));
words={ [fJ,fInv,fKp], [fKp,fInv,fJ], ...
    [fJ,fInv,fAk,fInv,fB], [fKp,fInv,fAj,fInv,fB], ...
    [fB,fInv,fAj,fInv,fKp], [fB,fInv,fAk,fInv,fJ], ...
    [fB,fInv,fAk,fInv,fAj,fInv,fB], ...
    [fB,fInv,fAj,fInv,fAk,fInv,fB] };
value=sym(0);
for n=1:numel(words)
    if ell==0, value=value+word(words{n}); else, value=value+wordq(words{n}); end
end
end

function f=factor(base,derivative)
f=struct('base',base,'derivative',derivative);
end

function value=select(vector,index)
if index==0, value=sym(0); else, value=vector(index); end
end

function value=word(factors)
value=sym(1);
for n=1:numel(factors), value=value*factors(n).base; end
end

function value=wordq(factors)
value=sym(0);
for active=1:numel(factors)
    term=factors(active).derivative;
    for n=1:numel(factors)
        if n~=active, term=term*factors(n).base; end
    end
    value=value+term;
end
end

function values=upper_double(exactValues)
values=zeros(size(exactValues));
for n=1:numel(exactValues), values(n)=exact_to_upper_double(exactValues(n)); end
end

function [value,index]=exact_min(values)
value=values(1); index=1;
for n=2:numel(values)
    if isAlways(values(n)<value), value=values(n); index=n; end
end
end

function out=failed_result(z,set)
out=struct('status','FAILED','bound',[],'linear_part',[],'remainder_part',[], ...
    'hessian_upper_bound',[],'parameter_box',set,'state',z,'runtime',0, ...
    'proof_method','','failure_reason','Uninitialized certificate.');
end
