function results=run_stage3a_tests()
%RUN_STAGE3A_TESTS Execute Stage-3A-specific MATLAB tests.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src'))); addpath(fullfile(root,'config'));
results=runtests(fullfile(root,'tests','test_hessian_certification.m'), ...
    'OutputDetail','Detailed');
assertSuccess(results);
end
