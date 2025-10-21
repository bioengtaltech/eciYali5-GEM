verbose = false;
saveModels = true;

% Get model adapter
adapterLocation = fullfile(pwd,'eciYali5GEMadapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

% Get pooled model path
eciYali5path = fullfile(params.path,'model','eciYali5-GEM.yml');

% Load pooled model
ecModel = readYAMLmodel(eciYali5path);

fluxData = loadFluxData();
lipidNchainData = loadLipidNchainData();
protData = loadProtData([3, 3, 3, 3, 3, 3]);

% Calculate GANnonPol for later
GAMnonPol = calculateGAM(ecModel);

for i = 1:length(fluxData.conds)
    % Update the biomass equation for each condition
    modelName = ['ec',fluxData.conds{i}];
    fprintf('Working on: %s \n', modelName)
    ecModel_new = ecModel;

    %% Update models genetic background
    if startsWith(fluxData.conds{i}, {'CAR'})

        % Add carotenoids
        metsToAdd.mets = {'phytoe','lycop','caro','NAD','NADH'};
        metsToAdd.metNames = {'phytoene','lycopene','beta-carotene','NAD','NADH'};
        metsToAdd.compartments = {'lp','lp','lp','lp','lp'};
        metsToAdd.unconstrained = [0, 0, 0, 0, 0];
        ecModel_new = addMets(ecModel_new,metsToAdd,true,'m');

        % Add heterologous genes
        genesToAdd.genes = {'CarRP','CarB'};
        genesToAdd.geneShortNames = {'CarRP','CarB'};
        ecModel_new = addGenesRaven(ecModel_new,genesToAdd);

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

        ecModel_new = addRxns(ecModel_new, rxnsToAdd, 3);

        [ecModel_new, addedRxns] = addExchangeRxns(ecModel_new,'out','caro');

        % Add enzymetic reactions
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

        newEnzymes.enzymes = {'Q9UUQ6','A0A168PH23'};
        newEnzymes.genes = {'CarRP', 'CarB'};
        newEnzymes.mw = [69841, 65633];

        ecModel_new = addNewRxnsToEC(ecModel_new, newRxns, newEnzymes, ModelAdapter);

        % Constrain the kcat
        % Calculated from specific activity in BRENDA (EC 2.5.1.32)
        ecModel_new = setKcatForReactions(ecModel_new,'PSY',0.0780);
        % Calculated from specific activity in BRENDA (EC 1.3.99.29)
        ecModel_new = setKcatForReactions(ecModel_new,'CRTI',4.3755e-04);
        % Calculated from specific activity in BRENDA (EC 5.5.1.19)
        ecModel_new = setKcatForReactions(ecModel_new,'LYCOPC',0.4144);

        ecModel_new = applyKcatConstraints(ecModel_new);

    elseif startsWith(fluxData.conds{i}, {'OBE'})

        % Gene deletions in OBE
        deletionsOBE = {'O74934',... % POX1
            'O74935',... % POX2
            'O74936',... % POX3
            'F2Z627',... % POX4
            'F2Z630',... % POX5
            'Q6C6T0',... % POX6
            'Q6C282'};   % TGL4
        for j = 1:length(deletionsOBE)
            rxns2knock = getReactionsFromEnzyme(ecModel,deletionsOBE{i});
            ecModel_new = setParam(ecModel,'eq',rxns2knock,0);
        end
    end


    %% Update Biomass Equation
    % Change protein content in the biomass equation
    ecModel_new = scaleBioMassYali4(ecModel_new, 'protein', fluxData.Ptot(i));

    % Change lipid content in the biomass equation
    ecModel_new = scaleBioMassYali4(ecModel_new, 'lipid', lipidNchainData.Ltot(i));

    % Adjust fatty acid distribution
    ecModel_new = updateAcylPool(ecModel_new, lipidNchainData, i, ModelAdapter);

    % Balance out mass with carbohydrate content
    [X,~,C,~,~,~,~] = sumBioMassYali4(ecModel_new, false);

    delta = X - 1;  % difference to balance
    fC = (C - delta)/C;
    ecModel_new = rescalePseudoReaction(ecModel_new, 'carbohydrate', fC);

    % Recalculate GAM
    [~,~,ecModel_new] = calculateGAM(ecModel_new, GAMnonPol,true);

    % Rename model
    ecModel_new.id = [modelName,'_pooled'];

    % Save pooled model
    if saveModels == true
        saveEcModel(ecModel_new,[ecModel_new.id,'.yml']);
    end

    % Constrain model with proteomics data
    ecModel_new = fillEnzConcs(ecModel_new, protData, i);
    ecModel_new = constrainEnzConcs(ecModel_new);

    % Update protein pool
    f = calculateFfactor(ecModel_new,protData);
    ecModel_new = setProtPoolSize(ecModel_new,fluxData.Ptot(i),f);

    feasability = -1;

    while feasability < 0
        % Load flux data
        ecModel_new = constrainFluxData(ecModel_new, fluxData, i,'max','loose');

        sol = solveLP(ecModel_new); % To observe if growth was reached.
        fprintf('Growth rate that is reached: %f /hour.\n', abs(sol.f))

        % This strain needs flexibilizeEnzConcs or it will not reach the
        % growth rate
        if startsWith(fluxData.conds{i},{'PARe','CARe','CARn','OBEe'})
            ecModel_new = setParam(ecModel_new,'lb','prot_pool_exchange',-1000);
        end
        fprintf('Flexibilizing protein bounds...\n')
        [ecModel_new, flexProt] = flexibilizeEnzConcs(ecModel_new, fluxData.grRate(i), 10, 5, ModelAdapter, verbose);
        fprintf('Flexibilization finished.\n')

        % New NGAM calc strategy
        ecModel_new = setParam(ecModel_new,'obj','xMAINTENANCE',1);
        sol = solveLP(ecModel_new,1);
        ecModel_new = setParam(ecModel_new,'lb','xMAINTENANCE',0.99*sol.f);
        ecModel_new = setParam(ecModel_new,'obj','xBIOMASS',1);

        % Flexibilize protein concentrations
        fprintf('Flexibilizing protein bounds...\n')
        [ecModel_new, flexProt] = flexibilizeEnzConcs(ecModel_new, fluxData.grRate(i), 10, 5, ModelAdapter, verbose);
        fprintf('Flexibilization finished.\n')
        % Sometime NGAM is too high and stops growth rate from being
        % achievable. This should fix that.
        sol = solveLP(ecModel_new);
        if sol.f < 0.975*fluxData.grRate(i)
            ecModel_new = setParam(ecModel_new,'lb','xMAINTENANCE',0);
            fprintf('Flexibilizing protein bounds...\n')
            [ecModel_new, flexProt] = flexibilizeEnzConcs(ecModel_new, fluxData.grRate(i), 10, 5, ModelAdapter, verbose);
            fprintf('Flexibilization finished.\n')
        end

        % Check if model works with tight contraints
        % Load flux data
        ecModel_new = constrainFluxData(ecModel_new, fluxData, i,'max',10);
        ecModel_new = setParam(ecModel_new,'ub','1992',0);
        ecModel_new = setParam(ecModel_new,'lb','1992',-1000);
        if startsWith(fluxData.conds{i},{'CARe','CARn'})
            ecModel_new = setParam(ecModel_new,'lb','1672',0);
        end
        sol = solveLP(ecModel_new);


        if sol.stat == 1

            ecModel_new = setParam(ecModel_new,'lb','xMAINTENANCE',0);
            sol = solveLP(ecModel_new);
            maxGrowth = sol.f;
            ecModel_new = setParam(ecModel_new,'lb','xBIOMASS',maxGrowth);
            ecModel_new = setParam(ecModel_new,'obj','xMAINTENANCE',1);
            sol = solveLP(ecModel_new,1);
            ecModel_new = setParam(ecModel_new,'lb','xMAINTENANCE',sol.f);
            break
        end
        feasability = sol.stat;
    end

    % Reset model protein pool and Rename model
    if startsWith(fluxData.conds{i},{'PARe','CARe','CARn','OBEe'})
        ecModel_new = setParam(ecModel_new,'lb','xBIOMASS',maxGrowth);
        ecModel_new = setParam(ecModel_new,'obj','prot_pool_exchange',1);
        sol = solveLP(ecModel_new);
        ecModel_new = setParam(ecModel_new,'lb','prot_pool_exchange',sol.f);
        ecModel_new = setParam(ecModel_new,'lb','xBIOMASS',0);
        ecModel_new = setParam(ecModel_new,'obj','xBIOMASS',1);
    end

    ecModel_new.id = [modelName,'_prot'];

    % Save proteome model
    if saveModels == true
        saveEcModel(ecModel_new,[ecModel_new.id,'.yml']);
    end
end