function tests=test_uncertainty_bounds
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
addpath(fullfile(root,'config'));
testCase.TestData.cfg=stage2_config();
end

function testParameterBoxes(testCase)
cfg=testCase.TestData.cfg;
for radius=cfg.relative_radii
    s=parameter_set(cfg.geometry,radius);
    verifyTrue(testCase,all(s.lower>0));
    verifyFalse(testCase,s.is_hardware_confidence_set);
    verifyFalse(testCase,s.includes_joint_physical_constraints);
end
end

function testCertifiedBoundContainsSamples(testCase)
cfg=testCase.TestData.cfg; rng(cfg.seed,'twister');
for radius=cfg.relative_radii
    s=parameter_set(cfg.geometry,radius);
    for stateIndex=[1,3,6]
        z=cfg.states(stateIndex);
        certificate=certified_reference_bound(z,s,cfg.geometry);
        verifyEqual(testCase,certificate.status,'CERTIFIED');
        F0=reference_torque_map(z,s.nominal,cfg.geometry);
        points=[s.lower,s.upper,s.lower+(s.upper-s.lower).*rand(6,12)];
        for k=1:size(points,2)
            actual=abs(reference_torque_map(z,points(:,k),cfg.geometry)-F0);
            verifyLessThanOrEqual(testCase,actual,certificate.channel_bound);
        end
    end
end
end

function testCertifiedComponents(testCase)
cfg=testCase.TestData.cfg; s=parameter_set(cfg.geometry,0.10); z=cfg.states(4);
c=certified_reference_bound(z,s,cfg.geometry);
verifyGreaterThanOrEqual(testCase,c.inertia_component,0);
verifyGreaterThanOrEqual(testCase,c.coriolis_component,0);
verifyEqual(testCase,c.channel_bound(1),c.inertia_component+c.coriolis_component, ...
    'RelTol',1e-14);
end

function testSchurAndInverseIdentities(testCase)
cfg=testCase.TestData.cfg; z=cfg.states(3); s=parameter_set(cfg.geometry,0.10);
theta=s.lower+0.37*(s.upper-s.lower);
p0=inertial_parameters_to_model(s.nominal,cfg.geometry);
p1=inertial_parameters_to_model(theta,cfg.geometry);
Q0=[0;0;0;z.q]; Qdot=[0;0;0;z.qdot];
b0=partition_dynamics(Q0,Qdot,p0); b1=partition_dynamics(Q0,Qdot,p1);
verifyEqual(testCase,b1.Mr,b1.D-b1.Bc.'*(b1.A\b1.Bc),'AbsTol',1e-13);
lhs=inv(b1.A)-inv(b0.A);
rhs=-inv(b0.A)*(b1.A-b0.A)*inv(b1.A);
verifyLessThan(testCase,norm(lhs-rhs,2),1e-13);
certificate=schur_error_bound(s,cfg.geometry);
verifyLessThanOrEqual(testCase,norm(b1.Mr-b0.Mr,2),certificate.epsilon_M);
end

function testCoriolisEnvelopeAtFixedState(testCase)
cfg=testCase.TestData.cfg; z=cfg.states(6); s=parameter_set(cfg.geometry,0.20);
theta=s.lower+0.63*(s.upper-s.lower);
p0=inertial_parameters_to_model(s.nominal,cfg.geometry);
p1=inertial_parameters_to_model(theta,cfg.geometry);
C0=reduced_coriolis_matrix(z.q,z.qdot,p0);
C1=reduced_coriolis_matrix(z.q,z.qdot,p1);
schur=schur_error_bound(s,cfg.geometry);
cor=coriolis_error_bound(z.qdot,schur);
verifyLessThanOrEqual(testCase,norm(C1-C0,2),cor.epsilon_C);
end

function testTaylorNeverSelfCertifies(testCase)
cfg=testCase.TestData.cfg; s=parameter_set(cfg.geometry,0.05);
t=taylor_remainder_bound(cfg.states(1),s,cfg.geometry,2,cfg.seed);
verifyEqual(testCase,t.status,'CANDIDATE');
verifyNotEqual(testCase,t.status,'CERTIFIED');
end
