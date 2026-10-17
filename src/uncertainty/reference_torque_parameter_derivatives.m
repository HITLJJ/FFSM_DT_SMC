function out=reference_torque_parameter_derivatives(z,theta,geometry)
%REFERENCE_TORQUE_PARAMETER_DERIVATIVES Analytic Schur parameter derivatives.
% Returns the 2x6 Jacobian and 6x6x2 Hessian of F [N*m] with respect to
% theta=[kg,kg,kg,kg*m^2,kg*m^2,kg*m^2]. Schur parameter derivatives are
% analytic; their joint-coordinate derivatives use complex step. This
% routine supports independent validation but is not itself a certificate.

theta=theta(:); q=z.q(:); h=1e-20;
[Mr,Mrj,Mrjk]=mass_derivatives(q,theta,geometry);
dMrj=zeros(2,2,2,6); dMrjk=zeros(2,2,2,6,6);
for ell=1:2
    qc=complex(q); qc(ell)=qc(ell)+1i*h;
    [~,Jc,Hc]=mass_derivatives(qc,theta,geometry);
    dMrj(:,:,ell,:)=imag(Jc)/h;
    dMrjk(:,:,ell,:,:)=imag(Hc)/h;
end

Cj=zeros(2,2,6); Cjk=zeros(2,2,6,6);
for a=1:2
    for b=1:2
        for ell=1:2
            for j=1:6
                Cj(a,b,j)=Cj(a,b,j)+0.5*(dMrj(a,b,ell,j)+ ...
                    dMrj(a,ell,b,j)-dMrj(b,ell,a,j))*z.qdot(ell);
                for k=1:6
                    Cjk(a,b,j,k)=Cjk(a,b,j,k)+0.5*( ...
                        dMrjk(a,b,ell,j,k)+dMrjk(a,ell,b,j,k)- ...
                        dMrjk(b,ell,a,j,k))*z.qdot(ell);
                end
            end
        end
    end
end

J=zeros(2,6); H=zeros(6,6,2);
for j=1:6
    J(:,j)=Mrj(:,:,j)*z.qr_ddot(:)+Cj(:,:,j)*z.qr_dot(:);
    for k=1:6
        H(j,k,:)=reshape(Mrjk(:,:,j,k)*z.qr_ddot(:)+ ...
            Cjk(:,:,j,k)*z.qr_dot(:),1,1,2);
    end
end
out.Mr=Mr; out.Mr_j=Mrj; out.Mr_jk=Mrjk;
out.C_j=Cj; out.C_jk=Cjk;
out.jacobian=J; out.hessian=H;
out.status='NUMERICAL_CROSSCHECK';
out.method='Analytic Schur parameter differentiation plus complex-step q differentiation.';
end

function [Mr,Mrj,Mrjk]=mass_derivatives(q,theta,geometry)
b=affine_mass_block_derivatives(q,theta,geometry);
A=b.A; B=b.B; D=b.D; K=A\eye(3,'like',A);
Mr=D-B.'*K*B;
Mrj=zeros(2,2,6,'like',A); Mrjk=zeros(2,2,6,6,'like',A);
for j=1:6
    Aj=b.Aj(:,:,j); Bj=b.Bj(:,:,j); Dj=b.Dj(:,:,j);
    Mrj(:,:,j)=Dj-Bj.'*K*B-B.'*K*Bj+B.'*K*Aj*K*B;
    for k=1:6
        Ak=b.Aj(:,:,k); Bk=b.Bj(:,:,k);
        S=Bj.'*K*Bk+Bk.'*K*Bj ...
            -Bj.'*K*Ak*K*B-Bk.'*K*Aj*K*B ...
            -B.'*K*Aj*K*Bk-B.'*K*Ak*K*Bj ...
            +B.'*K*Ak*K*Aj*K*B+B.'*K*Aj*K*Ak*K*B;
        Mrjk(:,:,j,k)=-S;
    end
end
end
