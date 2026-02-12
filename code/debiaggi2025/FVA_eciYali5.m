% STEP 73-75 Perform (ec)FVA
% Perform FVA on a conventional GEM, ecModel, and ecModel plus proteomics
% integration, all under similar exchange flux constraints. First make sure
% that the correct models are loaded.
model = loadConventionalGEM();
ecModel = loadEcModel('ecYeastGEM.yml');
ecModelProt = loadEcModel('ecYeastGEM_stage4.yml');

% As protein model can maximum reach 0.088, also set this as constrain for
% all models.
fluxData.grRate(1) = 0.0880;

% Apply same constraints on exchange fluxes
model = constrainFluxData(model,fluxData,1,'max','loose');
ecModel = constrainFluxData(ecModel,fluxData,1,'max','loose');
ecModelProt = constrainFluxData(ecModelProt,fluxData,1,'max','loose');

solveLP(model)
solveLP(ecModel)
solveLP(ecModelProt)

% Prepare output structure.
minFlux = zeros(numel(model.rxns),3);
maxFlux = minFlux;

% Run ecFVA for each model.
[minFlux(:,1), maxFlux(:,1)] = ecFVA(model, model);
[minFlux(:,2), maxFlux(:,2)] = ecFVA(ecModel, model);
[minFlux(:,3), maxFlux(:,3)] = ecFVA(ecModelProt, model);

% Plot ecFVA results and store in output/.
plotEcFVA(minFlux, maxFlux);
if savePlot == true
    saveas(gca, fullfile(params.path,'output','ecFVA.pdf'))
end