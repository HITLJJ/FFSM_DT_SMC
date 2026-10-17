function tests=test_reference_torque_map
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
addpath(fullfile(root,'config'));
testCase.TestData.p=default_parameters();
testCase.TestData.z=stage2_config().states(3);
end

function testDefinitionAndDimensions(testCase)
p=testCase.TestData.p; z=testCase.TestData.z;
[F,d]=reference_torque_map(z,nominal_inertial_parameters(p),p);
verifySize(testCase,F,[2,1]);
verifyEqual(testCase,F,d.Mr*z.qr_ddot+d.C*z.qr_dot,'AbsTol',1e-13);
end

function testNotPhysicalVelocitySubstitution(testCase)
p=testCase.TestData.p; z=testCase.TestData.z;
[F,d]=reference_torque_map(z,nominal_inertial_parameters(p),p);
wrong=d.Mr*z.qr_ddot+d.C*z.qdot;
verifyEqual(testCase,wrong-F,d.C*(z.qdot-z.qr_dot),'AbsTol',1e-13);
verifyGreaterThan(testCase,norm(wrong-F),1e-8);
end

function testModelParameterSeparation(testCase)
p=testCase.TestData.p; z=testCase.TestData.z;
t0=nominal_inertial_parameters(p); t1=t0; t1(2)=1.1*t1(2);
F0=reference_torque_map(z,t0,p); F1=reference_torque_map(z,t1,p);
verifyGreaterThan(testCase,norm(F1-F0),1e-8);
verifyEqual(testCase,nominal_inertial_parameters(p),t0);
end

function testInvalidParameterRejected(testCase)
p=testCase.TestData.p; z=testCase.TestData.z;
t=nominal_inertial_parameters(p); t(1)=-1;
verifyError(testCase,@()reference_torque_map(z,t,p),'FFSM:InvalidInertialParameters');
end
