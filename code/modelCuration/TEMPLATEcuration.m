%% TEMPLATE FOR CURATIONS APPLIED AFTER A RELEASE
% Copy this file to a versioned script, for example v1_0_1.m, and add all
% curations that should be applied on top of that released model.
%
% Recommended execution:
%   addpath(genpath('code'))
%   model = runCurationScript('v1_0_1', 'ModelName', 'eciYali5-GEM', 'Save', true);
%
% Curation scripts are run with these workspace variables available:
%   model      - model structure loaded by loadEciYaliModel
%   repoRoot   - repository root path
%   scriptName - curation script name

if ~exist('model', 'var')
    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(genpath(fullfile(repoRoot, 'code')));
    model = loadEciYaliModel('eciYali5-GEM', 'RepositoryRoot', repoRoot);
    scriptName = mfilename;
end

%% Brief description of curation (Issue/PR #xxx)
% Explain the biological or technical reason for the curation and cite
% supporting evidence. Keep any supporting tables in data/modelCuration/.

% Example:
% model = changeGeneAssoc(model, 'rxn_id', 'YALI0XXXXXg');

%% Validation
validateEciYaliModel(model);
