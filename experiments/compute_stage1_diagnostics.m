function d = compute_stage1_diagnostics(result, p, torqueFun)
%COMPUTE_STAGE1_DIAGNOSTICS Compute all stage-1 acceptance signals.

t = result.t;
x = result.x;
n = numel(t);
d.P = zeros(n,2);
d.LO = zeros(n,1);
d.T = zeros(n,1);
d.power = zeros(n,1);
d.qdd_full = zeros(n,2);
d.qdd_reduced = zeros(n,2);
d.acceleration_error = zeros(n,1);
d.base_velocity_reconstruction_error = zeros(n,1);
d.min_eig_H = zeros(n,1);
d.min_eig_Mr = zeros(n,1);
d.symmetry_H = zeros(n,1);
d.tau = zeros(n,2);

for k = 1:n
    Q = x(k,1:5).';
    Qdot = x(k,6:10).';
    tau = torqueFun(t(k),x(k,:).');
    me = momentum_energy(Q,Qdot,p);
    H = full_mass_matrix(Q,p);
    Mr = reduced_mass_matrix(Q(4:5),p);
    qddFull = full_acceleration(Q,Qdot,tau,p);
    qddReduced = reduced_acceleration(Q(4:5),Qdot(4:5),tau,p);
    reconstructed = zero_momentum_base_velocity(Q,Qdot(4:5),p);

    d.P(k,:) = me.P.';
    d.LO(k) = me.LO;
    d.T(k) = me.T;
    d.tau(k,:) = tau.';
    d.power(k) = Qdot(4:5).'*tau;
    d.qdd_full(k,:) = qddFull(4:5).';
    d.qdd_reduced(k,:) = qddReduced.';
    d.acceleration_error(k) = norm(qddFull(4:5)-qddReduced);
    d.base_velocity_reconstruction_error(k) = norm(Qdot(1:3)-reconstructed);
    d.min_eig_H(k) = min(eig((H+H.')/2));
    d.min_eig_Mr(k) = min(eig((Mr+Mr.')/2));
    d.symmetry_H(k) = norm(H-H.',2);
end

d.delta_P = d.P-d.P(1,:);
d.delta_LO = d.LO-d.LO(1);
if isfield(result,'actuator_work')
    d.work = result.actuator_work;
else
    % Backward-compatible fallback for previously saved result structures.
    d.work = cumtrapz(t,d.power);
end
d.work_energy_error = (d.T-d.T(1))-d.work;
d.max_linear_momentum_drift = max(vecnorm(d.delta_P.',2,1));
d.max_angular_momentum_drift = max(abs(d.delta_LO));
d.max_energy_drift = max(abs(d.T-d.T(1)));
d.max_work_energy_error = max(abs(d.work_energy_error));
d.max_acceleration_error = max(d.acceleration_error);
d.max_reconstruction_error = max(d.base_velocity_reconstruction_error);
d.max_symmetry_error = max(d.symmetry_H);
d.minimum_eigenvalue_H = min(d.min_eig_H);
d.minimum_eigenvalue_Mr = min(d.min_eig_Mr);
end
