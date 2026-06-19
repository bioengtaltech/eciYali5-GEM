function versions = readModelVersions(modelName, varargin)
% readModelVersions Read model-level release provenance metadata.
%
%   versions = readModelVersions()
%   versions = readModelVersions('eciYali5-GEM')

if nargin < 1 || isempty(modelName)
    modelName = '';
end

parser = inputParser;
addOptional(parser, 'modelName', modelName, @(x) ischar(x) || isstring(x));
addParameter(parser, 'RepositoryRoot', findRepoRoot(), @(x) ischar(x) || isstring(x));
parse(parser, modelName, varargin{:});

repoRoot = char(parser.Results.RepositoryRoot);
versionsFile = fullfile(repoRoot, 'model', 'versions.tsv');
if ~isfile(versionsFile)
    error('readModelVersions:MissingFile', 'Model provenance file not found: %s', versionsFile);
end

options = detectImportOptions(versionsFile, ...
    'FileType', 'text', ...
    'Delimiter', '\t', ...
    'TextType', 'string');
versions = readtable(versionsFile, options);

modelName = string(parser.Results.modelName);
if strlength(modelName) > 0
    versions = versions(strcmp(versions.model, modelName), :);
end
end

function repoRoot = findRepoRoot()
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(fileparts(thisFile)));
end
