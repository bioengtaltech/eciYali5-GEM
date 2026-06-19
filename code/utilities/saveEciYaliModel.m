function outputFiles = saveEciYaliModel(model, modelName, varargin)
% saveEciYaliModel Export an eciYali5-GEM model with repository defaults.
%
%   saveEciYaliModel(model)
%   saveEciYaliModel(model, 'iYali5', 'IncludeBinary', true)
%
% By default this writes XML/YML/TXT files only. Set IncludeBinary to true
% only for main/release workflows that intentionally need MAT/XLSX outputs.

if nargin < 2 || isempty(modelName)
    modelName = 'eciYali5-GEM';
end

parser = inputParser;
addParameter(parser, 'RepositoryRoot', findRepoRoot(), @(x) ischar(x) || isstring(x));
addParameter(parser, 'OutputDir', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Formats', {'xml', 'yml', 'txt'}, @(x) iscell(x) || isstring(x) || ischar(x));
addParameter(parser, 'IncludeBinary', false, @(x) islogical(x) || isnumeric(x));
parse(parser, varargin{:});

repoRoot = char(parser.Results.RepositoryRoot);
formats = normalizeFormats(parser.Results.Formats);
if parser.Results.IncludeBinary
    formats = unique([formats, {'xlsx', 'mat'}], 'stable');
end

[defaultDir, fileStem] = resolveModelLocation(repoRoot, char(modelName));
outputDir = char(parser.Results.OutputDir);
if isempty(outputDir)
    outputDir = defaultDir;
end

if exist('exportForGit', 'file') ~= 2
    error('saveEciYaliModel:MissingDependency', ...
        'exportForGit is required. Add RAVEN/GECKO utilities to the MATLAB path.');
end

if ~isfolder(outputDir)
    mkdir(outputDir);
end

exportForGit(model, fileStem, outputDir, formats, false, false);

outputFiles = cell(size(formats));
for i = 1:numel(formats)
    outputFiles{i} = fullfile(outputDir, [fileStem '.' formats{i}]);
end
end

function formats = normalizeFormats(formats)
if ischar(formats) || isstring(formats)
    formats = cellstr(formats);
end
formats = cellfun(@(x) lower(erase(char(x), '.')), formats, 'UniformOutput', false);
end

function [modelDir, fileStem] = resolveModelLocation(repoRoot, modelName)
switch lower(modelName)
    case {'eciyali5-gem', 'eciyali5'}
        modelDir = fullfile(repoRoot, 'model');
        fileStem = 'eciYali5-GEM';
    case {'iyali5', 'iyali5-gem'}
        modelDir = fullfile(repoRoot, 'model', 'iYali5-GEM');
        fileStem = 'iYali5';
    otherwise
        error('saveEciYaliModel:UnknownModel', ...
            'Unknown model "%s". Use "eciYali5-GEM" or "iYali5".', modelName);
end
end

function repoRoot = findRepoRoot()
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(fileparts(thisFile)));
end
