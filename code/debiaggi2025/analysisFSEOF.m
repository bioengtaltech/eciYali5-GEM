% Get model adapter
adapterLocation = fullfile(pwd,'eciYali5GEMadapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();
filePath = fullfile(params.path,'code','debiaggi2025','output','FSEOF');

% Load flux data
fluxData = loadFluxData();

% Conditionals
saveModels = true;

% Load original model
modelOG = loadConventionalGEM();

for i = 1:length(fluxData.conds)
    modelName = ['ec',fluxData.conds{i}];
    fprintf('Working on: %s \n', modelName)
    model = modelOG;

    % Load ecModel
    ecModel = readYAMLmodel(fullfile(params.path,'model','contextSpecific_eciYali5-GEM',[modelName,'_pooled.yml']));

    % Constrain models around experimental condition
    ecModel = setParam(ecModel,'eq','1808',fluxData.exchFluxes(i,1));
    model = setParam(model,'eq','1808',fluxData.exchFluxes(i,1));


    if startsWith(fluxData.conds{i}, {'CAR'})
        % ecModel FSEOF
        ecfseof = FSEOFalt(ecModel,'ecGEM','EXC_OUT_caro','1808',16,true,true,fluxData.conds{i},filePath);

        % Generate regular model for carotenoid production
        metsToAdd.mets = {'phytoe','lycop','caro','NAD','NADH'};
        metsToAdd.metNames = {'phytoene','lycopene','beta-carotene','NAD','NADH'};
        metsToAdd.compartments = {'lp','lp','lp','lp','lp'};
        metsToAdd.unconstrained = [0, 0, 0, 0, 0];
        model = addMets(model,metsToAdd,true,'m');

        % Add heterologous genes
        genesToAdd.genes = {'CarRP','CarB'};
        genesToAdd.geneShortNames = {'CarRP','CarB'};
        model = addGenesRaven(model,genesToAdd);

        % Add non enzymetic reactions
        rxnsToAdd.rxns = {'NADtlp',...
            'NADHtlp',...
            'GGDPtlp'};
        rxnsToAdd.equations = {'NAD[c] <=> NAD[lp]',...
            'NADH[c] <=> NADH[lp]',...
            'geranylgeranyl diphosphate[c] <=> geranylgeranyl diphosphate[lp]'};
        rxnsToAdd.rxnNames = {'NAD transport, cytoplasm-lipid particle',...
            'NADH transport, cytoplasm-lipid particle',...
            'geranylgeranyl diphosphate transport, cytoplasm-lipid particle'};
        rxnsToAdd.subSystems = {'transport, lipid particle',...
                                'transport, lipid particle',...
                                'transport, lipid particle'}';

        model = addRxns(model, rxnsToAdd, 3);

        [model, addedRxns] = addExchangeRxns(model,'out','caro');
        model.subSystems(end) = model.subSystems(end-1);

        % Add enzymatic reactions
        newRxns.rxns = {'PSY',...
            'CRTI',...
            'LYCOPC'};
        newRxns.rxnNames = {'phytoene synthase',...
            'phytoene dehydrogenase (lumped)',...
            'lycopene cyclase (lumped);'};
        newRxns.equations = {'geranylgeranyl diphosphate[lp] => phytoene[lp] + 2 diphosphate[lp]',...
            'phytoene[lp] + 4 NAD[lp] => lycopene[lp] + 4 NADH[lp]',...
            'lycopene[lp] => beta-carotene[lp]'};
        newRxns.grRules = {'CarRP',...
            'CarB',...
            'CarRP'};
        newRxns.subSystems = {'Carotenoid biosynthesis',...
                              'Carotenoid biosynthesis',...
                              'Carotenoid biosynthesis'}';

        model = addRxns(model, newRxns, 3);

        % Save model
        if saveModels == true
            exportModel(model,fullfile(params.path,'model','iYali5-GEM','iYali5_carotene.xml'));
        end

        % FSEOF for regular models
        fseof = FSEOFalt(model,'GEM','EXC_OUT_caro','1808',16,true,true,fluxData.conds{i},filePath);

    elseif startsWith(fluxData.conds{i}, {'PAR'}) || startsWith(fluxData.conds{i}, {'OBE'})
        ecfseof = FSEOFalt(ecModel,'ecGEM','EXC_OUT_m1640','1808',16,true,true,fluxData.conds{i},filePath);

        % FSEOF for regular models
        fseof = FSEOFalt(model,'GEM','EXC_OUT_m1640','1808',16,true,true,fluxData.conds{i},filePath);
    end
end