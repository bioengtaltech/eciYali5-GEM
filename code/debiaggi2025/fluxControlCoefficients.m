clear
clc

% Fetch this script's directory
scriptDir = fileparts(mfilename('fullpath'));

% Define output folder
baseOutputDir = fullfile(scriptDir,'output');

subDirPath = fullfile(baseOutputDir,'figures');

if ~exist(subDirPath, 'dir')
    mkdir(subDirPath);
end

adapterLocation = fullfile(pwd,'eciYali5GEMadapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

% Define data folder
dataDirPath = fullfile(params.path,'data');

% Load acronyms
acronymsTable = readtable(fullfile(dataDirPath,'rxnAcronyms.tsv'), 'FileType','delimitedtext','Delimiter', '\t');
names_data = acronymsTable.rxnNames;

acronymDictionary = struct();
acronymDictionary.rxnNames = acronymsTable.rxnNames;
acronymDictionary.acronyms = acronymsTable.Acronym;

modelFCC.names = {'ecCARe_pooled.yml','ecOBEe_pooled.yml'};
modelFCC.obj = {'EXC_OUT_caro','EXC_OUT_m1640'};
modelFCC.target = {'carotenoid','lipid'};

for i = 1:length(modelFCC.names)
    % Load the ecModel
    ecModel = readYAMLmodel(fullfile(params.path,'model','contextSpecific_eciYali5-GEM',modelFCC.names{i}));

    % Get relevant rxn indexes
    poolIdx = strcmpi(ecModel.rxns, 'prot_pool_exchange');

    % Maximize target, store original ecModel
    ecModel = setParam(ecModel, 'obj', modelFCC.obj{i}, 1);
    ecModel_OG = ecModel;
    sol = solveLP(ecModel);

    % Save reference target
    targetRef = sol.f;

    % Set Qa
    Qa = 1.001;

    % Initialize FCCs
    FCCs = zeros(numel(ecModel.ec.kcat), 1);

    progressbar('Flux control coefficients calculation')
    for j = 1:numel(ecModel.ec.kcat)
        % Calculate new kcat
        ecModel.ec.kcat(j) = ecModel.ec.kcat(j) * Qa;

        % Apply kcat constraints and solve
        ecModel = applyKcatConstraints(ecModel);
        sol = solveLP(ecModel, 1);
        targetFCC = sol.f;

        % Calculate FCC
        kcat_j = ecModel_OG.ec.kcat(j);
        FCCs(j) = (targetFCC - targetRef) * kcat_j / (targetRef * (kcat_j * Qa - kcat_j));

        % Restore the original ecModel
        ecModel = ecModel_OG;
        progressbar(j / numel(ecModel.ec.kcat))
    end
    progressbar(1)

    % Create output variable
    result.FCCs = FCCs/sum(FCCs);
    result.kcat = ecModel_OG.ec.kcat;
    result.rxns = ecModel_OG.ec.rxns;

    % Sort the results in descending order based on FCCs
    [result.FCCs, sortOrder] = sort(result.FCCs,'descend');
    result.kcat = result.kcat(sortOrder);
    result.rxns = result.rxns(sortOrder);
    result.rxnEnzMat = ecModel.ec.rxnEnzMat(sortOrder, :);

    % Initialize rxnNames and enzymes as cell arrays
    result.rxnNames = cell(length(result.rxns), 1);
    result.enzymes = cell(length(result.rxns), 1);

    % Populate rxnNames and enzymes
    for j = 1:length(result.rxns)
        result.rxnNames(j,1) = ecModel.rxnNames(strcmp(ecModel.rxns,result.rxns(j)));
        if j < length(ecModel.ec.enzymes)
            result.enzymes{j} = ecModel.ec.enzymes(logical(result.rxnEnzMat(j,:))');
            if length(result.enzymes{j}) > 1
                result.enzymes{j} = strjoin(result.enzymes{j},', ');
            end
        end
    end

%% Make chart

% Sample data
FCCs = result.FCCs(1:10); % Top 10 FCCs
rxnNames = result.rxnNames(1:10); % Top 10 rxnNames
kcats = result.kcat(1:10);

acronyms4plot = cell(size(rxnNames));
for j = 1:length(rxnNames)
    current_name = rxnNames{j};
    idx = find(strcmp(acronymDictionary.rxnNames,current_name));
    acronyms4plot{j} = acronymDictionary.acronyms{idx};
end

% Colors
pleasantColors = [0.2 0.4 0.8; 0.8 0.2 0.4; 0.4 0.6 0.2; 0.6 0.2 0.6; 0.2 0.6 0.6;
                  0.8 0.4 0.2; 0.4 0.2 0.8; 0.6 0.6 0.2; 0.2 0.8 0.4; 0.4 0.4 0.4];

%% ---- Figure setup ----
figWidth_mm  = 80;                   % target figure width
figHeight_mm = figWidth_mm * 3/4;    % 4:3 aspect ratio
innerWidth_cm  = 6;                  % fixed inner axes width
innerHeight_cm = 4.5;                % fixed inner axes height

fig = figure('Units', 'centimeters');
ax = axes('Parent', fig, 'Units', 'centimeters');

% Fix data area size (inner plotting area)
ax.Position = [2 1.5 innerWidth_cm innerHeight_cm];
ax.ActivePositionProperty = 'position'; % ensures this position is respected

%% ---- Plot ----
h = barh(ax, FCCs);
set(ax, 'YDir', 'reverse', 'XScale', 'log');
yticks(1:length(acronyms4plot));
yticklabels(acronyms4plot);

xlabel('Flux Control Coefficients', 'FontName', 'Arial', 'FontSize', 10);
ylabel('Enzyme (k_{cat} s^{-1})', 'FontName', 'Arial', 'FontSize', 10);

grid(ax, 'off');
set(ax, 'FontName', 'Arial', 'FontSize', 10);
set(ax, 'Box', 'on', 'LineWidth', 0.75, 'Color', 'none');
ax.XAxis.Color = 'black';
ax.YAxis.Color = 'black';
ax.XAxis.LineWidth = 0.75;
ax.YAxis.LineWidth = 0.75;
ax.XAxis.TickDirection = 'both';
ax.YAxis.TickDirection = 'both';

% Bar colors
for j = 1:numel(h)
    set(h(j), 'FaceColor', pleasantColors(j,:), 'EdgeColor', 'none');
end

% Text labels
for j = 1:length(FCCs)
    text(FCCs(j) * 1.05, j, sprintf('%.2e', FCCs(j)), ...
        'VerticalAlignment', 'middle', 'FontSize', 10, 'FontName', 'Arial');
end

%% ---- Adjust figure to prevent clipping ----
ti = ax.TightInset; % [left bottom right top]
outerWidth = innerWidth_cm + ti(1) + ti(3) + 2;
outerHeight = innerHeight_cm + ti(2) + ti(4) + 2;
set(fig, 'Position', [2 2 outerWidth outerHeight]);

%% ---- Save as SVG ----
figureName = ['FCC_', ecModel.id, '_', modelFCC.target{i}, '.svg'];
saveas(fig, fullfile(subDirPath, figureName), 'svg');

end