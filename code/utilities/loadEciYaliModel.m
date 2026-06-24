function [model, modelFile] = loadEciYaliModel(modelName, varargin)
% loadEciYaliModel Load an eciYali5-GEM repository model file.
%
%   model = loadEciYaliModel()
%   model = loadEciYaliModel('eciYali5-GEM')
%   model = loadEciYaliModel('iYali5')
%   model = loadEciYaliModel(..., 'Format', 'xml')
%
% This mirrors the Yeast-GEM loader pattern: YML/YAML files are loaded with
% readYAMLmodel, other text model formats use RAVEN importModel when
% available, and COBRA readCbModel is used only as a fallback. MAT files are
% loaded directly and must contain a model structure.

if nargin < 1 || isempty(modelName)
    modelName = 'eciYali5-GEM';
end

parser = inputParser;
addOptional(parser, 'modelName', 'eciYali5-GEM', @(x) ischar(x) || isstring(x));
addParameter(parser, 'RepositoryRoot', findRepoRoot(), @(x) ischar(x) || isstring(x));
addParameter(parser, 'Format', 'auto', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Path', '', @(x) ischar(x) || isstring(x));
parse(parser, modelName, varargin{:});

modelName = char(parser.Results.modelName);
repoRoot = char(parser.Results.RepositoryRoot);
requestedFormat = lower(char(parser.Results.Format));
explicitPath = char(parser.Results.Path);

if ~isempty(explicitPath)
    modelFile = explicitPath;
else
    [modelDir, fileStem] = resolveModelLocation(repoRoot, modelName);
    modelFile = resolveModelFile(modelDir, fileStem, requestedFormat);
end

if ~isfile(modelFile)
    error('loadEciYaliModel:MissingFile', 'Model file not found: %s', modelFile);
end

[~, ~, ext] = fileparts(modelFile);
switch lower(ext)
    case '.mat'
        model = loadMatModel(modelFile);
    case {'.yml', '.yaml'}
        if exist('readYAMLmodel', 'file') ~= 2
            error('loadEciYaliModel:MissingDependency', ...
                ['readYAMLmodel is required to load YAML model files. ' ...
                'Add RAVEN to the MATLAB path.']);
        end
        model = readYAMLmodel(modelFile);
    otherwise
        if ~(exist('ravenCobraWrapper.m', 'file') == 2)
            if exist('readCbModel.m', 'file') == 2
                warning(['RAVEN cannot be found. eciYali5-GEM will instead be loaded in ' ...
                    'COBRA format.\n\nRAVEN is recommended when curating eciYali5-GEM.'])
                model = readCbModel(modelFile);
            else
                error('loadEciYaliModel:MissingDependency', ...
                    'RAVEN cannot be found. See README.md for installation instructions.');
            end
        else
            model = importModel(modelFile);
        end
end
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
        error('loadEciYaliModel:UnknownModel', ...
            'Unknown model "%s". Use "eciYali5-GEM" or "iYali5".', modelName);
end
end

function modelFile = resolveModelFile(modelDir, fileStem, requestedFormat)
if strcmp(requestedFormat, 'auto')
    formats = {'yml', 'xml', 'txt', 'mat'};
else
    formats = {erase(requestedFormat, '.')};
end

modelFile = '';
for i = 1:numel(formats)
    candidate = fullfile(modelDir, [fileStem '.' formats{i}]);
    if isfile(candidate)
        modelFile = candidate;
        return
    end
end

error('loadEciYaliModel:MissingFile', ...
    'No %s model file found in %s.', strjoin(formats, '/'), modelDir);
end

function model = loadMatModel(modelFile)
loaded = load(modelFile);
if isfield(loaded, 'model')
    model = loaded.model;
    return
end

names = fieldnames(loaded);
for i = 1:numel(names)
    candidate = loaded.(names{i});
    if isstruct(candidate) && all(isfield(candidate, {'rxns', 'mets', 'S'}))
        model = candidate;
        return
    end
end

error('loadEciYaliModel:MissingModelStruct', ...
    'MAT file does not contain a recognizable model structure: %s', modelFile);
end

function repoRoot = findRepoRoot()
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(fileparts(thisFile)));
end
