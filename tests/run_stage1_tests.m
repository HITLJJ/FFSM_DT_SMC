function results = run_stage1_tests()
%RUN_STAGE1_TESTS Execute all automated stage-1 MATLAB tests.
% Entry point: addpath(genpath(pwd)); results = run_stage1_tests;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
addpath(fullfile(root,'config'));
addpath(fullfile(root,'experiments'));
results = runtests(fileparts(mfilename('fullpath')), ...
    'IncludeSubfolders',false,'OutputDetail','Detailed');
assertSuccess(results);
end
