function ref=closed_loop_reference(t)
%CLOSED_LOOP_REFERENCE Smooth feasible two-joint reference in rad and seconds.
q0=[0.40;-0.70]; amplitude=[0.22;-0.18]; omega=[0.75;0.60];
ref.q=q0+amplitude.*(1-cos(omega*t));
ref.qdot=amplitude.*omega.*sin(omega*t);
ref.qddot=amplitude.*omega.^2.*cos(omega*t);
end
