function [GAMnonPol, GAM, model] = calculateGAM(model, GAMnonPol, setGAM, modelAdapter)
    % calculateGAM
    %   Calculate either GAM or GAMnonPol based on the input model
    %   and optionally set GAM in the model
    %
    % Input:
    %   model           - The metabolic model
    %   GAMnonPol       - Optional, only if you want to calculate GAMnonPol
    %   setGAM          - Optional, default is false
    %   modelAdapter    - Optional, model adapter (default is loaded with ModelAdapterManager)
    %
    % Output: 
    %   GAMnonPol       - Calculated GAMnonPol
    %   GAM             - Calculated GAM
    %   model           - model with updated GAM
    %

    % Default values
    if nargin < 4 || isempty(modelAdapter)
        modelAdapter = ModelAdapterManager.getDefault();
        if isempty(modelAdapter)
            error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
        end
    end

    if nargin < 3
        setGAM = false;
    end

    if nargin < 2
        GAMnonPol = [];
    end

    % Get model parameters
    params = modelAdapter.getParameters();

    % Calculate GAMpol
    [~, P, C, R, D, ~, ~] = sumBioMassYali5(model, false);
    GAMpol = P * 37.7 + C * 12.8 + R * 26.0 + D * 26.0;

    % Calculate GAMnonPol if not provided
    if isempty(GAMnonPol)
        bioRxnIdx = getIndexes(model, params.bioRxn, 'rxns');
        ADP = find(strcmp(model.metNames, 'ADP'));
        GAEC = max(full(model.S(ADP, bioRxnIdx)));
        GAMnonPol = GAEC - GAMpol;
    end

    % Calculate GAM
    
    GAM = GAMpol + GAMnonPol;

    % Set GAM in the model if requested
    if setGAM
        bioRxnIdx = getIndexes(model, params.bioRxn, 'rxns');
        for i = 1:length(model.mets)
            S_ix = model.S(i, bioRxnIdx);
            isGAM = sum(strcmp({'ATP', 'ADP', 'H2O', 'H+_p+1', 'phosphate'}, model.metNames{i})) == 1; % those are metNames specific to the iYali_corr model
            if S_ix ~= 0 && isGAM
                model.S(i, bioRxnIdx) = sign(S_ix) * GAM;
            end
        end
    end
end