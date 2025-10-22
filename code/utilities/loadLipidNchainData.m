function lipidNchainData = loadLipidNchainData(lipidNchainFile, modelAdapter)
% loadLipidNchainData

if nargin < 2 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if nargin < 1 || isempty(lipidNchainFile)
    lipidNchainFile = fullfile(params.path,'data','lipidNchainData.tsv');
end

% Load lipids composition and fatty acid molar composition data:

fid = fopen(lipidNchainFile);

lipidNchainComp = textscan(fid, '%s %f32 %f32 %f32 %f32 %f32 %f32', 'Delimiter', '\t', 'HeaderLines', 1);

lipidNchainData.conds = lipidNchainComp{1};
lipidNchainData.Ltot = lipidNchainComp{2};

lipidNchainData.chainConds = [];

for j = 3:length(lipidNchainComp)
    lipidNchainData.chainConds = [lipidNchainData.chainConds, lipidNchainComp{j}];
end
lipidNchainData.chainConds = double(lipidNchainData.chainConds);

fclose(fid);
end