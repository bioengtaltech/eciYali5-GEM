clear
clc

% Fetch this script's directory
scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);

% Save outputs in R_analysis folder
Routputs_path = fullfile(scriptDir,'R_Analysis','data','processed');

% Define output folder
baseOutputDir = fullfile(scriptDir,'output');

% Create results subdirectories
subDirs = 'escher_jsons';
subDirPath = fullfile(baseOutputDir,subDirs);

if ~exist(subDirPath, 'dir')
    mkdir(subDirPath);
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

for i = 1:length(fluxData.conds)

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

    % Get index of C-source rxn
    idx = getIndexes(ecModel,params.c_source,'rxns');

    % Compute absolute value of the element in the specified row
    abs_value = abs(solutions(idx, :));

    % Divide each column by the absolute value
    solutionsNorm = solutions ./ abs_value;

    % Store solutions in a field named after the model name
    allSolutions.(modelName) = solutionsNorm;
    allSol4Plot.(modelName) = mapRxnsToConv(ecModel,model,solutionsNorm);
    allSol4PCA.(modelName) = full(mean(solutionsNorm,2));

    % Data treatment
    fluxMean = full(mean(solutionsNorm,2));
    fluxSD = full(std(solutionsNorm,0,2));
    logVec = abs(fluxMean) > abs(fluxSD);
    fluxMeanSD = logVec.*fluxMean;

    % Export data
    fluxesForEscher(ecModel.rxns,fluxMean,[modelName,'_prot_FBA_SD.json'],subDirPath);
    fluxesForEscher(model.rxns,mapRxnsToConv(ecModel,model,fluxMean),[modelName,'_FBA_SD.json'],subDirPath);

end
%% Flexible Mean Calculation and Comparison
% List of condition pairs
conditionPairs = {
    'ecPARe', 'ecPARn';
    'ecCARe', 'ecCARn';
    'ecOBEe', 'ecOBEn';
    'ecPARe', 'ecCARe';
    'ecPARe', 'ecOBEe';
    'ecPARn', 'ecCARn';
    'ecPARn', 'ecOBEn';
    };

% Threshold for solving precision
threshold = 1e-8;

% Loop through each pair of conditions
for pairIdx = 1:size(conditionPairs, 1)
    condition1 = conditionPairs{pairIdx, 1};
    condition2 = conditionPairs{pairIdx, 2};

    % Calculate means of the fluxes for the selected conditions
    mean_fluxes_condition1 = mean(full(allSol4Plot.(condition1)), 2);
    mean_fluxes_condition1(abs(mean_fluxes_condition1) < threshold) = 0;

    mean_fluxes_condition2 = mean(full(allSol4Plot.(condition2)), 2);
    mean_fluxes_condition2(abs(mean_fluxes_condition2) < threshold) = 0;

    % Calculate fold changes (or differences in means)
    log2FC = log2(abs(mean_fluxes_condition2) ./ abs(mean_fluxes_condition1));

    % Perform t-tests across the two conditions
    [~, p_values] = ttest2(full(allSol4Plot.(condition2)), full(allSol4Plot.(condition1)), 'Dim', 2);

    % Adjust p-values using Benjamini-Hochberg FDR
    adjusted_p_values = mafdr(p_values, 'BHFDR', true);

    % Generate CSV File with Results
    % Define the output filename (unique for each pair)
    outputFilename = sprintf('flux_analysis_results_%s_vs_%s.csv',condition2, condition1);
    outputFilename = fullfile(Routputs_path, outputFilename);

    % Add 'r_' prefix to reaction identifiers
    rxnIds = strcat('r_', model.rxns);

    % Create a table with the required columns
    resultsTable = table(rxnIds, log2FC, adjusted_p_values, ...
        'VariableNames', {'Entry', 'log2FC', 'adj_p_value'});

    % Write the table to a CSV file with a semicolon as the separator
    writetable(resultsTable, outputFilename, 'Delimiter', ';');
    disp(['Results saved to ', outputFilename]);
end

%% Datasets creation for R analysis

% Create annotations and rxnGoTerms tables (done once for the entire model)
% Define the annotations output filename
annotationsTableName = fullfile(projectRoot,'debiaggi2025','R_Analysis','data','raw','subSystemsAnnotations.csv');

% Create a table with the required columns
annotationsTable = table(model.rxnNames, model.subSystems, ...
    'VariableNames', {'Entry','subSystem'});

% Write the table to a CSV file
writetable(annotationsTable, annotationsTableName, 'Delimiter', ';');
disp(['Annotations saved to ', annotationsTableName]);

% Create a table with prefixed reaction IDs and gene-reaction rules
rxnGoTermsTable = table(rxnIds, model.grRules, ...
    'VariableNames', {'rxns', 'genes'});

% Export the table to a CSV file (semicolon-delimited)
rxnGoTermsFilename = fullfile(projectRoot,'debiaggi2025','R_Analysis','data','raw','rxnGoTermsTable.csv');
writetable(rxnGoTermsTable, rxnGoTermsFilename, 'Delimiter', ';');
disp(['Reaction GO Terms saved to ', rxnGoTermsFilename]);