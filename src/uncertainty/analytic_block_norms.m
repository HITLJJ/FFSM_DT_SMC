function coeff = analytic_block_norms(geometry, exactMode)
%ANALYTIC_BLOCK_NORMS Global rational norm bounds for H block coefficients.
% Bounds are uniform over all q. They follow from |sin|,|cos|<=1,
% ||X'Y||_2<=||X||_F||Y||_F, and ab<=(a^2+b^2)/2. No state samples enter.
% Rows correspond to [m0,m1,m2,J0,J1,J2]; derivative columns to q1,q2.

if nargin<2, exactMode=false; end
if exactMode
    cv=@(x) sym(x,'f');
else
    cv=@(x) x;
end
a0=cv(geometry.a0); c1=cv(geometry.c1);
l1=cv(geometry.l1); c2=cv(geometry.c2);

jb1=2+(a0+c1)^2; jq1=c1^2;
jb2=2+(a0+l1+c2)^2;
jq2=(l1+c2)^2+c2^2;

coeff.A=[2;jb1;jb2;1;1;1];
coeff.B=[0;(jb1+jq1)/2;(jb2+jq2)/2;0;1;2];
coeff.D=[0;jq1;jq2;0;1;2];

z=cv(0);
coeff.Aq=repmat(z,6,2); coeff.Bq=repmat(z,6,2); coeff.Dq=repmat(z,6,2);
% Link 1, q1 derivative.
djb=c1^2; djq=c1^2;
coeff.Aq(2,1)=jb1+djb;
coeff.Bq(2,1)=(djb+jq1+jb1+djq)/2;
coeff.Dq(2,1)=jq1+djq;
% Link 2, q1 derivative.
djb=(l1+c2)^2; djq=jq2;
coeff.Aq(3,1)=jb2+djb;
coeff.Bq(3,1)=(djb+jq2+jb2+djq)/2;
coeff.Dq(3,1)=jq2+djq;
% Link 2, q2 derivative.
djb=c2^2; djq=2*c2^2;
coeff.Aq(3,2)=jb2+djb;
coeff.Bq(3,2)=(djb+jq2+jb2+djq)/2;
coeff.Dq(3,2)=jq2+djq;
end
