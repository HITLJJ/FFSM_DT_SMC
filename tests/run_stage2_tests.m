function results=run_stage2_tests()
%RUN_STAGE2_TESTS Execute Stage-2-specific MATLAB tests.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
addpath(fullfile(root,'config'));
files={fullfile(root,'tests','test_reference_torque_map.m'), ...
    fullfile(root,'tests','test_uncertainty_bounds.m')};
results=runtests(files,'OutputDetail','Detailed');
assertSuccess(results);
end
