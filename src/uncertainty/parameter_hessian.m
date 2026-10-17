function H = parameter_hessian(z, theta, geometry)
%PARAMETER_HESSIAN Numerical parameter Hessians of both F channels.
% H is 6x6x2. Values are exploratory and must not be used as global certified
% Hessian envelopes without an independent validated enclosure.

theta=theta(:); n=6; H=zeros(n,n,2);
h=8e-4*max(abs(theta),1);
f0=reference_torque_map(z,theta,geometry);
for i=1:n
    di=zeros(n,1); di(i)=h(i);
    fp=reference_torque_map(z,theta+di,geometry);
    fm=reference_torque_map(z,theta-di,geometry);
    H(i,i,:)=reshape((fp-2*f0+fm)/h(i)^2,1,1,2);
    for j=i+1:n
        dj=zeros(n,1); dj(j)=h(j);
        val=(reference_torque_map(z,theta+di+dj,geometry) ...
            -reference_torque_map(z,theta+di-dj,geometry) ...
            -reference_torque_map(z,theta-di+dj,geometry) ...
            +reference_torque_map(z,theta-di-dj,geometry))/(4*h(i)*h(j));
        H(i,j,:)=reshape(val,1,1,2); H(j,i,:)=reshape(val,1,1,2);
    end
end
end
