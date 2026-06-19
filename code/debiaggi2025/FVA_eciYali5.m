clear
clc

% Fetch this script's directory
scriptDir = fileparts(mfilename('fullpath'));

% Get model adapter
adapterLocation = fullfile(pwd,'eciYali5GEMadapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

% Define output folder
outputDir = fullfile(scriptDir,'output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Create images subdirectory
imagesDir = fullfile(outputDir, 'figures');
if ~exist(imagesDir, 'dir')
    mkdir(imagesDir);
end

% Define paths for models
specificEcGEMpath = fullfile(params.path,'model','contextSpecific_eciYali5-GEM');
baseEcGEMpath = fullfile(params.path,'model','eciYali5-GEM.yml');

% Load flux data and conventional model
fluxData = loadFluxData();
model = loadConventionalGEM();
lipidNchainData = loadLipidNchainData();

% Load base ecModel
ecModel_base = readYAMLmodel(baseEcGEMpath);

% Calculate GAMnonPol for later
GAMnonPol = calculateGAM(model);

for i = 1:length(fluxData.conds)
    
    condName = fluxData.conds{i};
    modelName = ['ec', condName];
    fprintf('Performing FVA for: %s \n', modelName)
    
    % Load specific proteome model
    fileName = [modelName,'_prot.yml'];
    protModelPath = fullfile(specificEcGEMpath,fileName);
    
    if ~exist(protModelPath, 'file')
        warning('Model %s not found. Skipping...', fileName);
        continue;
    end
    
    ecModel_prot = readYAMLmodel(protModelPath);
    
    % Apply constraints
    % 1. Conventional GEM
    model_cond = model;
    model_cond = scaleBioMassYali5(model_cond, 'protein', fluxData.Ptot(i));
    model_cond = scaleBioMassYali5(model_cond, 'lipid', lipidNchainData.Ltot(i));
    model_cond = updateAcylPool(model_cond, lipidNchainData, i, ModelAdapter);
    [X,~,C,~,~,~,~] = sumBioMassYali5(model_cond, false);
    delta = X - 1;
    fC = (C - delta)/C;
    model_cond = rescalePseudoReaction(model_cond, 'carbohydrate', fC);
    [~,~,model_cond] = calculateGAM(model_cond, GAMnonPol, true);
    model_cond = constrainFluxData(model_cond, fluxData, i, 'max', 'loose');
    
    % 2. Base ecModel
    ecModel_base_cond = constrainFluxData(ecModel_base, fluxData, i, 'max', 'loose');
    
    % 3. Proteome ecModel
    ecModel_prot_cond = constrainFluxData(ecModel_prot, fluxData, i, 'max', 'loose');

    
    if startsWith(fluxData.conds{i},{'CARn'})
            model_cond = setParam(model_cond,'lb','1992',-1000);
            ecModel_base_cond = setParam(ecModel_base_cond,'lb','1992',-1000);
            ecModel_prot_cond = setParam(ecModel_prot_cond,'lb','1992',-1000);
    end
    % % Equalize growth rates
    % sol1 = solveLP(model_cond);
    % sol2 = solveLP(ecModel_base_cond);
    % sol3 = solveLP(ecModel_prot_cond);
    % 
    % minGR = min([abs(sol1.f), abs(sol2.f), abs(sol3.f)]);
    % 
    % model_cond = setParam(model_cond, 'lb', params.bioRxn, minGR);
    % ecModel_base_cond = setParam(ecModel_base_cond, 'lb', params.bioRxn, minGR);
    % ecModel_prot_cond = setParam(ecModel_prot_cond, 'lb', params.bioRxn, minGR);
    
    % Prepare output structure
    minFlux = zeros(numel(model.rxns),3);
    maxFlux = minFlux;
    
    % Run ecFVA for each model
    fprintf('Running FVA on conventional GEM...\n');
    [minFlux(:,1), maxFlux(:,1)] = ecFVA(model_cond, model_cond);
    
    fprintf('Running FVA on base ecModel...\n');
    [minFlux(:,2), maxFlux(:,2)] = ecFVA(ecModel_base_cond, model_cond);
    
    fprintf('Running FVA on proteome-constrained ecModel...\n');
    [minFlux(:,3), maxFlux(:,3)] = ecFVA(ecModel_prot_cond, model_cond);
    
    % Plot ecFVA results
    plotEcFVA(minFlux, maxFlux);
    title('');
    
    % Save plot
    saveas(gca, fullfile(imagesDir, ['ecFVA_', condName, '.svg']));
    close(gcf);
    
end