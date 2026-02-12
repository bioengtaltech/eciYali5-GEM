%% Checking if ENO is a limiting enzyme:
clear
clc

% Fetch this script's directory
scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);


% Define output folder
baseOutputDir = fullfile(scriptDir,'output');

% Create enzymeUsage subdirectory
enzymeUsageDir = fullfile(baseOutputDir, 'enzymeUsage');
if ~exist(enzymeUsageDir, 'dir')
    mkdir(enzymeUsageDir);
end

% Get model adapter and load flux data
adapterLocation = fullfile(pwd,'eciYali5GEMadapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

% Set specific models path
specificEcGEMpath = fullfile(params.path,'model','contextSpecific_eciYali5-GEM');

fluxData = loadFluxData();
model = loadConventionalGEM();
allSolutions = struct();  % Initialize allSolutions as a structure
allUsageReports = struct(); % Initialize allUsageReports as a structure
allRobustness = struct(); % Initialize allRobustness as a structure

% Run analyzeEcGEMs pipeline, but just for ecCARe and ecCARn
for i = 3:4

    modelName = ['ec',fluxData.conds{i}];
    fprintf('Working on: %s \n', modelName)

    fileName = [modelName,'_prot.yml'];

    % Load model: already constrained with 10% variance around chemostat fluxes

    ecModel = readYAMLmodel(fullfile(specificEcGEMpath,fileName));

    % Set bounds for biomass too
    ecModel = setParam(ecModel,'obj','xBIOMASS',1);
    sol = solveLP(ecModel);
    ecModel = setParam(ecModel,'var', 'xBIOMASS', sol.f, 10);

    % Minimize all protein usage

    % Minimize total enzyme usage
    ecModel = setParam(ecModel,'obj','prot_pool_exchange',1);
    sol = solveLP(ecModel);
    ecModel = setParam(ecModel,'var', 'prot_pool_exchange', sol.f, 10);

    if startsWith(fluxData.conds{i}, {'CAR'})
        ecModel = setParam(ecModel,'obj','EXC_OUT_caro',1);
        sol = solveLP(ecModel);
        ecModel = setParam(ecModel,'var', 'EXC_OUT_caro', sol.f, 10);
    end

    % Get good reactions
    [~, goodRxns] = randomSampling(ecModel,1,true,true,true);

    % Run random sampling
    solutions = randomSampling(ecModel,5000,true,true,true,goodRxns);

    % Store solutions in a field named after the model name
    allSolutions.(modelName) = solutions;

    % Inspect enzyme usage
    fprintf('Reporting enzyme usage for %s\n', modelName);
    usageData = enzymeUsage(ecModel, full(mean(solutions,2)));
    allUsageData.(modelName) = usageData;
    usageReport = reportEnzymeUsage(ecModel,usageData);
    allUsageReports.(modelName) = usageReport;
    fprintf('Top absolute enzyme usage for %s:\n', modelName);
    disp(usageReport.topAbsUsage);

    % Save usageReport tables
    fields = fieldnames(usageReport);
    for f = 1:length(fields)
        fieldName = fields{f};
        if istable(usageReport.(fieldName))
            writetable(usageReport.(fieldName), fullfile(enzymeUsageDir, [modelName, '_', fieldName, '.csv']));
        end
    end

    % Robustness analysis
    ecModel = setParam(ecModel,'lb','EXC_OUT_caro',0);
    ecModel = setParam(ecModel,'ub','EXC_OUT_caro',1000);
    if i == 4
        % To avoid infeasible solutions in runRobustnessAnalysis() for CARn
        ecModel = setParam(ecModel,'ub','usage_prot_Q6C1F3',-0.4577);
        [controlFlux, objFlux] = runRobustnessAnalysis(ecModel, 'usage_prot_Q6C1F3', 20, 'EXC_OUT_caro', true);
    end
    allRobustness.(modelName).controlFlux = controlFlux;
    allRobustness.(modelName).objFlux = objFlux;
end