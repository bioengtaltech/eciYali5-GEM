function [model, modelFile, worktreePath] = getEarlierModelVersion(refName, modelName, varargin)
% getEarlierModelVersion Load a model from a git ref without switching branches.
%
%   model = getEarlierModelVersion('v1.0.1')
%   [model, modelFile] = getEarlierModelVersion('main', 'iYali5')
%
% A temporary detached git worktree is created and removed automatically
% unless KeepWorktree is true.

if nargin < 2 || isempty(modelName)
    modelName = 'eciYali5-GEM';
end

parser = inputParser;
addRequired(parser, 'refName', @(x) ischar(x) || isstring(x));
addOptional(parser, 'modelName', modelName, @(x) ischar(x) || isstring(x));
addParameter(parser, 'RepositoryRoot', findRepoRoot(), @(x) ischar(x) || isstring(x));
addParameter(parser, 'Format', 'auto', @(x) ischar(x) || isstring(x));
addParameter(parser, 'KeepWorktree', false, @(x) islogical(x) || isnumeric(x));
parse(parser, refName, modelName, varargin{:});

repoRoot = char(parser.Results.RepositoryRoot);
worktreePath = tempname;

addCommand = sprintf('git -C "%s" worktree add --detach --quiet "%s" "%s"', ...
    repoRoot, worktreePath, char(parser.Results.refName));
[status, output] = system(addCommand);
if status ~= 0
    error('getEarlierModelVersion:GitWorktreeFailed', ...
        'Could not create worktree for %s:%s%s', char(parser.Results.refName), newline, output);
end

cleanupObj = [];
if ~parser.Results.KeepWorktree
    cleanupObj = onCleanup(@() removeWorktree(repoRoot, worktreePath)); %#ok<NASGU>
end

[model, modelFile] = loadEciYaliModel(char(parser.Results.modelName), ...
    'RepositoryRoot', worktreePath, ...
    'Format', char(parser.Results.Format));
end

function removeWorktree(repoRoot, worktreePath)
if isfolder(worktreePath)
    command = sprintf('git -C "%s" worktree remove --force "%s"', repoRoot, worktreePath);
    system(command);
end
end

function repoRoot = findRepoRoot()
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(fileparts(thisFile)));
end
