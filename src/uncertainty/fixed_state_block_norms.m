function coeff=fixed_state_block_norms(q,geometry)
%FIXED_STATE_BLOCK_NORMS Certified fixed-q block-coefficient norm bounds.
% Trigonometric values are enclosed with exact-rational Taylor polynomials
% and Lagrange remainders. Matrix interval products use exact endpoints.
% Returned spectral-norm bounds are rational upper bounds on Frobenius norms.

assert(license('test','Symbolic_Toolbox'),'FFSM:SymbolicUnavailable', ...
    'Symbolic Math Toolbox is required.');
q=sym(q(:),'f'); assert(numel(q)==2,'FFSM:Dimension','q must be 2x1.');
a0=sym(geometry.a0,'f'); l1=sym(geometry.l1,'f');
c1=sym(geometry.c1,'f'); c2=sym(geometry.c2,'f');
S=iv_exact([0,-1;1,0]); I=iv_exact(eye(2)); Z21=iv_exact(zeros(2,1));

[cQ1,sQ1]=trig_enclosure(q(1));
[cQ12,sQ12]=trig_enclosure(q(1)+q(2));
d0=iv_exact([a0;0]); d1c=iv_scale(iv_stack(cQ1,sQ1),c1);
d1=iv_scale(iv_stack(cQ1,sQ1),l1);
d2c=iv_scale(iv_stack(cQ12,sQ12),c2);

J0=iv_exact([eye(2),zeros(2,3)]);
J1=iv_hcat(I,iv_mul(S,iv_add(d0,d1c)),iv_mul(S,d1c),Z21);
J2=iv_hcat(I,iv_mul(S,iv_add(iv_add(d0,d1),d2c)), ...
    iv_mul(S,iv_add(d1,d2c)),iv_mul(S,d2c));

Jq=cell(3,2); Jq{1,1}=iv_exact(zeros(2,5)); Jq{1,2}=Jq{1,1};
Jq{2,1}=iv_hcat(iv_exact(zeros(2,2)),iv_neg(d1c),iv_neg(d1c),Z21);
Jq{2,2}=iv_exact(zeros(2,5));
Jq{3,1}=iv_hcat(iv_exact(zeros(2,2)),iv_neg(iv_add(d1,d2c)), ...
    iv_neg(iv_add(d1,d2c)),iv_neg(d2c));
Jq{3,2}=iv_hcat(iv_exact(zeros(2,2)),iv_neg(d2c),iv_neg(d2c),iv_neg(d2c));
J={J0,J1,J2};

Hw=cell(3,1);
w0=[0,0,1,0,0]; w1=[0,0,1,1,0]; w2=[0,0,1,1,1];
Hw{1}=iv_exact(w0.'*w0); Hw{2}=iv_exact(w1.'*w1); Hw{3}=iv_exact(w2.'*w2);
H=cell(6,1); Hq=cell(6,2);
for body=1:3
    H{body}=iv_mul(iv_transpose(J{body}),J{body});
    H{body+3}=Hw{body};
    for ell=1:2
        Hq{body,ell}=iv_add(iv_mul(iv_transpose(Jq{body,ell}),J{body}), ...
            iv_mul(iv_transpose(J{body}),Jq{body,ell}));
        Hq{body+3,ell}=iv_exact(zeros(5));
    end
end

coeff.A=sym(zeros(6,1)); coeff.B=sym(zeros(6,1)); coeff.D=sym(zeros(6,1));
coeff.Aq=sym(zeros(6,2)); coeff.Bq=sym(zeros(6,2)); coeff.Dq=sym(zeros(6,2));
for j=1:6
    coeff.A(j)=iv_norm(H{j},1:3,1:3);
    coeff.B(j)=iv_norm(H{j},1:3,4:5);
    coeff.D(j)=iv_norm(H{j},4:5,4:5);
    for ell=1:2
        coeff.Aq(j,ell)=iv_norm(Hq{j,ell},1:3,1:3);
        coeff.Bq(j,ell)=iv_norm(Hq{j,ell},1:3,4:5);
        coeff.Dq(j,ell)=iv_norm(Hq{j,ell},4:5,4:5);
    end
end
coeff.proof_method=['Exact-rational trigonometric Taylor enclosure and ', ...
    'exact endpoint matrix interval arithmetic at fixed q.'];
coeff.trig_terms=18;
end

function [c,s]=trig_enclosure(x)
terms=18; sx=sym(0); cx=sym(0);
for n=0:terms-1
    sx=sx+(-1)^n*x^(2*n+1)/factorial(sym(2*n+1));
    cx=cx+(-1)^n*x^(2*n)/factorial(sym(2*n));
end
% Lagrange remainder after degrees 2*terms-1 (sin) and 2*terms-2 (cos).
rs=abs(x)^(2*terms)/factorial(sym(2*terms));
rc=abs(x)^(2*terms-1)/factorial(sym(2*terms-1));
s=iv(sx-rs,sx+rs); c=iv(cx-rc,cx+rc);
end

function a=iv(lo,hi), a=struct('lo',lo,'hi',hi); end
function a=iv_exact(x)
if ~isa(x,'sym'), x=sym(x,'f'); end
a=iv(x,x);
end
function a=iv_stack(varargin)
lo=cellfun(@(x)x.lo,varargin,'UniformOutput',false);
hi=cellfun(@(x)x.hi,varargin,'UniformOutput',false);
a=iv(vertcat(lo{:}),vertcat(hi{:}));
end
function a=iv_hcat(varargin)
lo=cellfun(@(x)x.lo,varargin,'UniformOutput',false);
hi=cellfun(@(x)x.hi,varargin,'UniformOutput',false);
a=iv(horzcat(lo{:}),horzcat(hi{:}));
end
function c=iv_add(a,b), c=iv(a.lo+b.lo,a.hi+b.hi); end
function c=iv_neg(a), c=iv(-a.hi,-a.lo); end
function c=iv_transpose(a), c=iv(a.lo.',a.hi.'); end
function c=iv_scale(a,s)
assert(isAlways(s>=0),'FFSM:IntervalScale','Geometry scale must be nonnegative.');
c=iv(a.lo*s,a.hi*s);
end

function c=iv_mul(a,b)
assert(size(a.lo,2)==size(b.lo,1),'FFSM:IntervalDimension','Invalid interval product.');
c=iv(sym(zeros(size(a.lo,1),size(b.lo,2))),sym(zeros(size(a.lo,1),size(b.lo,2))));
for row=1:size(c.lo,1)
    for col=1:size(c.lo,2)
        for inner=1:size(a.lo,2)
            products=[a.lo(row,inner)*b.lo(inner,col),a.lo(row,inner)*b.hi(inner,col), ...
                a.hi(row,inner)*b.lo(inner,col),a.hi(row,inner)*b.hi(inner,col)];
            c.lo(row,col)=c.lo(row,col)+min(products);
            c.hi(row,col)=c.hi(row,col)+max(products);
        end
    end
end
end

function upper=iv_norm(a,rows,cols)
mag=max(abs(a.lo(rows,cols)),abs(a.hi(rows,cols)));
radicand=sum(mag.^2,'all');
upper=sym(exact_to_upper_double(sqrt(radicand)),'f');
end
