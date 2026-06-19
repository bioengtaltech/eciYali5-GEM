function report = validateEciYaliModel(model, varargin)
% validateEciYaliModel Run lightweight structural and solver checks.

parser = inputParser;
addRequired(parser, 'model', @isstruct);
addParameter(parser, 'RunSolver', true, @(x) islogical(x) || isnumeric(x));
parse(parser, model, varargin{:});

report = struct();
report.ok = true;
report.errors = {};
report.warnings = {};

requiredFields = {'rxns', 'mets', 'S'};
for i = 1:numel(requiredFields)
    if ~isfield(model, requiredFields{i})
        report = addError(report, sprintf('Missing required model field: %s', requiredFields{i}));
    end
end

if report.ok
    [nMets, nRxns] = size(model.S);
    if numel(model.mets) ~= nMets
        report = addError(report, 'Number of metabolites does not match rows in S.');
    end
    if numel(model.rxns) ~= nRxns
        report = addError(report, 'Number of reactions does not match columns in S.');
    end
    if numel(unique(model.mets)) ~= numel(model.mets)
        report = addWarning(report, 'Model contains duplicate metabolite IDs.');
    end
    if numel(unique(model.rxns)) ~= numel(model.rxns)
        report = addWarning(report, 'Model contains duplicate reaction IDs.');
    end
end

if report.ok && parser.Results.RunSolver
    report = runSolverCheck(model, report);
end

if ~report.ok
    error('validateEciYaliModel:InvalidModel', '%s', strjoin(report.errors, newline));
end
end

function report = runSolverCheck(model, report)
if exist('solveLP', 'file') ~= 2
    report = addWarning(report, 'solveLP is not available; solver check skipped.');
    return
end
if ~isfield(model, 'c') && ~isfield(model, 'obj')
    report = addWarning(report, 'No objective field found; solver check skipped.');
    return
end

try
    sol = solveLP(model, 1);
    if ~isstruct(sol) || ~isfield(sol, 'stat')
        report = addWarning(report, 'solveLP returned an unexpected result structure.');
    elseif sol.stat <= 0
        report = addWarning(report, sprintf('solveLP did not report an optimal solution (stat=%g).', sol.stat));
    end
catch err
    report = addWarning(report, sprintf('solveLP check skipped after error: %s', err.message));
end
end

function report = addError(report, message)
report.ok = false;
report.errors{end + 1} = message;
end

function report = addWarning(report, message)
report.warnings{end + 1} = message;
end
