function results=run_closed_loop_tests()
%RUN_CLOSED_LOOP_TESTS Execute Plant-DT-SMC interface and regression tests.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src'))); addpath(fullfile(root,'config'));
addpath(fullfile(root,'experiments'));
files={fullfile(root,'tests','test_smc_closed_loop.m'), ...
    fullfile(root,'tests','test_dt_update.m'), ...
    fullfile(root,'tests','test_closed_loop_regression.m')};
results=runtests(files,'OutputDetail','Detailed'); assertSuccess(results);
end
