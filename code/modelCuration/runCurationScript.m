function model = runCurationScript(curationScript, varargin)
% runCurationScript Load a model, run a curation script, validate and save.
%
%   model = runCurationScript('v1_0_1')
%   model = runCurationScript('v1_0_1', 'Save', true)

parser = inputParser;
addRequired(parser, 'curationScript', @(x) ischar(x) || isstring(x));
addParameter(parser, 'ModelName', 'eciYali5-GEM', @(x) ischar(x) || isstring(x));
addParameter(parser, 'RepositoryRoot', findRepoRoot(), @(x) ischar(x) || isstring(x));
addParameter(parser, 'BaseRef', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Save', false, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'SaveFormats', {'xml', 'yml', 'txt'}, @(x) iscell(x) || isstring(x) || ischar(x));
addParameter(parser, 'IncludeBinary', false, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'Validate', true, @(x) islogical(x) || isnumeric(x));
parse(parser, curationScript, varargin{:});

repoRoot = char(parser.Results.RepositoryRoot);
addpath(genpath(fullfile(repoRoot, 'code')));

if strlength(string(parser.Results.BaseRef)) > 0
    model = getEarlierModelVersion(char(parser.Results.BaseRef), char(parser.Results.ModelName), ...
        'RepositoryRoot', repoRoot);
else
    model = loadEciYaliModel(char(parser.Results.ModelName), 'RepositoryRoot', repoRoot);
end

scriptPath = resolveScriptPath(repoRoot, char(parser.Results.curationScript));
[~, scriptName] = fileparts(scriptPath); %#ok<ASGLU>
run(scriptPath);

if parser.Results.Validate
    validateEciYaliModel(model);
    assertDevelopModelPolicy('RepositoryRoot', repoRoot);
end

if parser.Results.Save
    saveEciYaliModel(model, char(parser.Results.ModelName), ...
        'RepositoryRoot', repoRoot, ...
        'Formats', parser.Results.SaveFormats, ...
        'IncludeBinary', parser.Results.IncludeBinary);
end
end

function scriptPath = resolveScriptPath(repoRoot, curationScript)
[folder, name, ext] = fileparts(curationScript);
if isempty(ext)
    ext = '.m';
end
if isempty(folder)
    scriptPath = fullfile(repoRoot, 'code', 'modelCuration', [name ext]);
else
    scriptPath = fullfile(folder, [name ext]);
end
if ~isfile(scriptPath)
    error('runCurationScript:MissingScript', 'Curation script not found: %s', scriptPath);
end
end

function repoRoot = findRepoRoot()
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(fileparts(thisFile)));
end
