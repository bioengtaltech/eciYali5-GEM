function forbiddenFiles = assertDevelopModelPolicy(varargin)
% assertDevelopModelPolicy Enforce develop-branch model binary policy.
%
% Generated model MAT/XLSX files must not be committed on develop. Source
% model files under model/**/sourceModels are allowed.

parser = inputParser;
addParameter(parser, 'RepositoryRoot', findRepoRoot(), @(x) ischar(x) || isstring(x));
addParameter(parser, 'Enforce', [], @(x) islogical(x) || isnumeric(x) || isempty(x));
parse(parser, varargin{:});

repoRoot = char(parser.Results.RepositoryRoot);
branchName = currentBranch(repoRoot);
enforce = parser.Results.Enforce;
if isempty(enforce)
    enforce = strcmp(branchName, 'develop');
end

forbiddenFiles = findForbiddenModelBinaries(repoRoot);
if enforce && ~isempty(forbiddenFiles)
    error('assertDevelopModelPolicy:ForbiddenBinaries', ...
        'Generated model binaries are not allowed on develop:%s%s', ...
        newline, strjoin(forbiddenFiles, newline));
end
end

function forbiddenFiles = findForbiddenModelBinaries(repoRoot)
modelRoot = fullfile(repoRoot, 'model');
files = [dir(fullfile(modelRoot, '**', '*.mat')); dir(fullfile(modelRoot, '**', '*.xlsx'))];
forbiddenFiles = {};

for i = 1:numel(files)
    if files(i).isdir
        continue
    end
    fullPath = fullfile(files(i).folder, files(i).name);
    relPath = relativePath(repoRoot, fullPath);
    normalized = strrep(relPath, '\', '/');
    if contains(normalized, '/sourceModels/')
        continue
    end
    forbiddenFiles{end + 1} = relPath; %#ok<AGROW>
end
end

function branchName = currentBranch(repoRoot)
[status, output] = system(sprintf('git -C "%s" branch --show-current', repoRoot));
if status == 0
    branchName = strtrim(output);
else
    branchName = '';
end
end

function relPath = relativePath(repoRoot, fullPath)
repoRoot = char(java.io.File(repoRoot).getCanonicalPath());
fullPath = char(java.io.File(fullPath).getCanonicalPath());
prefix = [repoRoot filesep];
if startsWith(fullPath, prefix)
    relPath = extractAfter(fullPath, strlength(prefix));
else
    relPath = fullPath;
end
end

function repoRoot = findRepoRoot()
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(fileparts(thisFile)));
end
