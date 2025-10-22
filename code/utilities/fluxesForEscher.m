function reactionData = fluxesForEscher(reactionNames, fluxValues, fileName, filePath, modelAdapter)
% makeFluxesForEscher - Generate a .json file of flux values for Escher visualization
%
% Syntax:
%   reactionData = makeFluxesForEscher(reactionNames, fluxValues, fileName, filePath, modelAdapter)
%
% Description:
%   This function creates a .json file containing flux values for a set of
%   reactions, suitable for visualization in Escher maps.
%
% Input:
%   - reactionNames: Cell array of strings representing the names of biochemical reactions of a genome-scale model.
%   - fluxValues: Numeric array of flux values corresponding to the provided reaction names.
%   - fileName: (Optional) Name of the output file (default: 'escherFluxes.json').
%   - filePath: (Optional) Path to the directory where the file will be saved (default: 'output' within the model's path).
%   - modelAdapter: (Optional) Adapter for the biochemical model (default: default model adapter).
%
% Output:
%   - reactionData: Structure containing the reaction names and their corresponding flux values.
%
% Usage:
%   1. Call the function with the required inputs to generate the Escher-compatible .json file.
%   2. The output structure 'reactionData' provides a convenient representation of the data.
%
% Notes:
%   - If the modelAdapter is not provided, the function uses the default model adapter. If none is set, an error is raised.
%   - The output file is saved with the specified fileName and filePath, ensuring the directory exists or creating it if necessary.
%   - The generated .json file is suitable for Escher visualization tools.
%
%
% Author: Juliano Sabedotti De Biaggi
% Created: 08-Jan-2024
% Version: 1.0
% Revision Date: [Revision Date]

if nargin < 5 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefault();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.getParameters();

if nargin < 4 || isempty(filePath)
    filePath = fullfile(params.path, 'output');
end

if nargin < 3 || isempty(fileName)
    fileName = 'escherFluxes.json';  % Change the file extension to .txt
end

% Check if reactionNames and fluxValues are of the same size
if length(reactionNames) ~= length(fluxValues)
    error('Error: The number of reaction names must be equal to the number of flux values.');
end

% Filter out reaction names with null (zero) flux values
nonNullIndices = fluxValues ~= 0;
reactionNames = reactionNames(nonNullIndices);
fluxValues = fluxValues(nonNullIndices);

% Ensure the directory exists, create it if needed
if ~exist(filePath, 'dir')
    mkdir(filePath);
end

% Full file path including directory
fullFilePath = fullfile(filePath, fileName);

% Open the file for writing
fid = fopen(fullFilePath, 'w');
if fid > 0
    % Populate the structure and write data to the file
    fprintf(fid, '{');
    for i = 1:length(reactionNames)
        fprintf(fid, '"%s": %f', reactionNames{i}, fluxValues(i));
        % Add a comma if it's not the last entry
        if i < length(reactionNames)
            fprintf(fid, ',');
        end        
    end
    fprintf(fid, '}');
    fclose(fid);
    disp(['Text file exported successfully to: ' fullFilePath]);
else
    error('Unable to open the file for writing.');
end
end