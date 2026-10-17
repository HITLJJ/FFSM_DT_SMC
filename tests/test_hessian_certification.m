function tests=test_hessian_certification
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src'))); addpath(fullfile(root,'config'));
cfg=stage3a_config(); z=cfg.states(3);
testCase.TestData.cfg=cfg; testCase.TestData.z=z;
testCase.TestData.coeff=fixed_state_block_norms(z.q,cfg.geometry);
end

function testAnalyticDerivativesAgainstIndependentDifferences(testCase)
cfg=testCase.TestData.cfg; z=testCase.TestData.z; theta=nominal_inertial_parameters(cfg.geometry);
analytic=reference_torque_parameter_derivatives(z,theta,cfg.geometry);
finiteJ=parameter_jacobian(z,theta,cfg.geometry);
finiteH=parameter_hessian(z,theta,cfg.geometry);
verifyLessThan(testCase,norm(analytic.jacobian-finiteJ,'fro'),1e-8);
verifyLessThan(testCase,norm(analytic.hessian-finiteH,'fro'),2e-6);
verifyLessThan(testCase,norm(analytic.hessian(:,:,1)-analytic.hessian(:,:,1).','fro'),1e-12);
verifyLessThan(testCase,norm(analytic.hessian(:,:,2)-analytic.hessian(:,:,2).','fro'),1e-12);
end

function testFixedStateBlockEnclosures(testCase)
cfg=testCase.TestData.cfg; z=testCase.TestData.z;
theta=nominal_inertial_parameters(cfg.geometry); c=testCase.TestData.coeff;
b=affine_mass_block_derivatives(z.q,theta,cfg.geometry);
for j=1:6
    verifyLessThanOrEqual(testCase,norm(b.Aj(:,:,j),2),double(c.A(j))+1e-12);
    verifyLessThanOrEqual(testCase,norm(b.Bj(:,:,j),2),double(c.B(j))+1e-12);
    verifyLessThanOrEqual(testCase,norm(b.Dj(:,:,j),2),double(c.D(j))+1e-12);
end
h=1e-6;
for ell=1:2
    d=zeros(2,1); d(ell)=h;
    bp=affine_mass_block_derivatives(z.q+d,theta,cfg.geometry);
    bm=affine_mass_block_derivatives(z.q-d,theta,cfg.geometry);
    for j=1:6
        verifyLessThanOrEqual(testCase,norm((bp.Aj(:,:,j)-bm.Aj(:,:,j))/(2*h),2), ...
            double(c.Aq(j,ell))+1e-7);
        verifyLessThanOrEqual(testCase,norm((bp.Bj(:,:,j)-bm.Bj(:,:,j))/(2*h),2), ...
            double(c.Bq(j,ell))+1e-7);
        verifyLessThanOrEqual(testCase,norm((bp.Dj(:,:,j)-bm.Dj(:,:,j))/(2*h),2), ...
            double(c.Dq(j,ell))+1e-7);
    end
end
end

function testCompleteCertificateFieldsAndCoverage(testCase)
cfg=testCase.TestData.cfg; z=testCase.TestData.z;
s=parameter_set(cfg.geometry,0.10);
c=certified_hessian_bound(z,s,cfg.geometry,testCase.TestData.coeff);
verifyEqual(testCase,c.status,'CERTIFIED');
required={'bound','linear_part','remainder_part','hessian_upper_bound', ...
    'parameter_box','state','runtime','proof_method','failure_reason'};
verifyTrue(testCase,all(isfield(c,required)));
rng(cfg.seed,'twister'); F0=reference_torque_map(z,s.nominal,cfg.geometry);
points=[s.lower,s.upper,s.lower+(s.upper-s.lower).*rand(6,10)];
for n=1:size(points,2)
    actual=abs(reference_torque_map(z,points(:,n),cfg.geometry)-F0);
    verifyLessThanOrEqual(testCase,actual,c.bound);
    derivative=reference_torque_parameter_derivatives(z,points(:,n),cfg.geometry);
    for channel=1:2
        verifyLessThanOrEqual(testCase,abs(derivative.hessian(:,:,channel)), ...
            c.hessian_upper_bound(:,:,channel));
    end
end
end

function testRemainderScalingAndSelection(testCase)
cfg=testCase.TestData.cfg; s=parameter_set(cfg.geometry,0.20);
c=certified_hessian_bound(testCase.TestData.z,s,cfg.geometry,testCase.TestData.coeff);
verifyEqual(testCase,c.status,'CERTIFIED');
verifyLessThanOrEqual(testCase,c.remainder_part,c.elementwise_remainder);
verifyLessThanOrEqual(testCase,c.remainder_part,c.raw_spectral_remainder);
verifyLessThanOrEqual(testCase,c.remainder_part,c.normalized_spectral_remainder);
verifyEqual(testCase,c.normalized_radius,0.2*ones(6,1),'AbsTol',1e-15);
end

function testNoSamplesInCertificate(testCase)
cfg=testCase.TestData.cfg; s=parameter_set(cfg.geometry,0.05);
c=certified_hessian_bound(testCase.TestData.z,s,cfg.geometry,testCase.TestData.coeff);
verifyFalse(testCase,c.uses_samples);
verifyEqual(testCase,c.failure_reason,'');
end
