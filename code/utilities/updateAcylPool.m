function ecModel = updateAcylPool(ecModel,lipidNchainData,cond,modelAdapter)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

% Load acyl chain metabolites
acylMetabPath = fullfile(params.path,'data','xPOOLmets.tsv');

fid = fopen(acylMetabPath);
data = textscan(fid, '%s %s %s %s %s %s %s', 'Delimiter', '\t', 'HeaderLines', 1);
fclose(fid);

chainData.rxns = data{1};

chainData.mets = [];

for i = 2:length(data)
    chainData.mets = [chainData.mets, data{i}];
end


% Get acyl Pool rxns
FArxns = ecModel.rxns(find(contains(ecModel.rxns,'xPOOL')));

for i = 1:length(FArxns)
    idx = find(strcmp(chainData.rxns,FArxns(i)));
    equations.mets = chainData.mets(idx,:);
    % Adjust the signal of the coefficients according to the direction of the
    % rxn
    if contains(FArxns(i),'_REV')
        equations.stoichCoeffs = -[-lipidNchainData.chainConds(cond,:),1];
    else
        equations.stoichCoeffs = [-lipidNchainData.chainConds(cond,:),1];
    end

    % Do the actual adjustment in the FA/Acyl-FA pools
    ecModel = changeRxns(ecModel,FArxns(i),equations,1);
end
end