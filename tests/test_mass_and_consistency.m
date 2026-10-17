function tests = test_mass_and_consistency
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
cfg = stage1_config();
testCase.TestData.p = cfg.parameters;
testCase.TestData.thresholds = cfg.thresholds;
rng(20261016,'twister');
end

function testMassSymmetryAndPositiveDefiniteness(testCase)
p = testCase.TestData.p; th = testCase.TestData.thresholds;
worstSymmetry = 0; minH = inf; minMr = inf;
for k = 1:100
    Q = [2*randn(2,1); pi*(2*rand(3,1)-1)];
    H = full_mass_matrix(Q,p);
    Mr = reduced_mass_matrix(Q(4:5),p);
    worstSymmetry = max(worstSymmetry,norm(H-H.',2));
    minH = min(minH,min(eig((H+H.')/2)));
    minMr = min(minMr,min(eig((Mr+Mr.')/2)));
end
verifyLessThanOrEqual(testCase,worstSymmetry,th.symmetry);
verifyGreaterThan(testCase,minH,th.positive_definite);
verifyGreaterThan(testCase,minMr,th.positive_definite);
end

function testIndependentKineticEnergy(testCase)
p = testCase.TestData.p;
for k = 1:25
    Q = randn(5,1); Qdot = randn(5,1);
    H = full_mass_matrix(Q,p);
    matrixEnergy = 0.5*Qdot.'*H*Qdot;
    directEnergy = momentum_energy(Q,Qdot,p).T;
    verifyLessThan(testCase,abs(matrixEnergy-directEnergy),1e-11);
end
end

function testNonlinearVectorAgainstFiniteDifferenceChristoffel(testCase)
p = testCase.TestData.p; th = testCase.TestData.thresholds;
fdStep = 2e-6;
for sample = 1:10
    Q = randn(5,1); Qdot = randn(5,1);
    dH = zeros(5,5,5);
    for k = 1:5
        delta = zeros(5,1); delta(k)=fdStep;
        dH(:,:,k) = (full_mass_matrix(Q+delta,p)-full_mass_matrix(Q-delta,p))/(2*fdStep);
    end
    hfd = zeros(5,1);
    for i = 1:5
        for j = 1:5
            for k = 1:5
                gamma = 0.5*(dH(i,j,k)+dH(i,k,j)-dH(j,k,i));
                hfd(i) = hfd(i)+gamma*Qdot(j)*Qdot(k);
            end
        end
    end
    h = full_nonlinear_vector(Q,Qdot,p);
    verifyLessThan(testCase,norm(h-hfd),th.finite_difference);
end
end

function testFullAndReducedAcceleration(testCase)
p = testCase.TestData.p; th = testCase.TestData.thresholds;
for k = 1:30
    Q = [randn(2,1); randn(3,1)];
    qdot = randn(2,1); tau = randn(2,1);
    bdot = zero_momentum_base_velocity(Q,qdot,p);
    Qdot = [bdot;qdot];
    qddFull = full_acceleration(Q,Qdot,tau,p);
    qddReduced = reduced_acceleration(Q(4:5),qdot,tau,p);
    verifyLessThan(testCase,norm(qddFull(4:5)-qddReduced),th.acceleration);
end
end

function testZeroMomentumInitialization(testCase)
p = testCase.TestData.p;
Q = [0.7;-0.4;0.3;0.6;-0.8]; qdot = [0.25;-0.18];
Qdot = [zero_momentum_base_velocity(Q,qdot,p);qdot];
me = momentum_energy(Q,Qdot,p);
verifyLessThan(testCase,norm(me.P),1e-12);
verifyLessThan(testCase,abs(me.LO),1e-12);
end

function testOptionalPayloadCountedOnce(testCase)
p = testCase.TestData.p; Q = randn(5,1);
H0 = full_mass_matrix(Q,p);
p.payload.enabled = true; p.payload.model = 'point';
p.payload.mass = 0.75; p.payload.J = 0;
H1 = full_mass_matrix(Q,p);
kin = planar_kinematics(Q,p);
Jp = kin.Jv{end};
expected = p.payload.mass*(Jp.'*Jp);
verifyLessThan(testCase,norm((H1-H0)-expected,2),1e-11);
end
