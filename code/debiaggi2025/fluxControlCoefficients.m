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
        acronyms4plot{j} = [acronymDictionary.acronyms{idx}];
    end

    % Define a set of pleasant colors
    pleasantColors = [0.2 0.4 0.8; 0.8 0.2 0.4; 0.4 0.6 0.2; 0.6 0.2 0.6; 0.2 0.6 0.6;
        0.8 0.4 0.2; 0.4 0.2 0.8; 0.6 0.6 0.2; 0.2 0.8 0.4; 0.4 0.4 0.4];

    % Create a horizontal bar chart
    h = barh(FCCs);

    % Customize the chart
    xlabel('Flux Control Coefficients');
    ylabel('Enzyme (k_{cat} s^{-1})');

    % Set custom labels for the y-axis ticks
    yticks(1:10); % Set the number of ticks to match the number of bars
    yticklabels(acronyms4plot); % Set the rxnNames as labels

    % Set the chart area size to 5 in. tall by 6.75 in. wide
    set(gcf, 'Position', [100, 100, 2 * 6.75 * 100, 5 * 100]); % [left, bottom, width, height]

    % Invert the y-axis to display the highest value at the top
    set(gca, 'YDir', 'reverse');

    % Set x-axis to logarithmic scale
    set(gca, 'XScale', 'log');

    % Remove grid lines
    grid off;

    % Set font size for labels and ticks
    set(gca, 'FontSize', 14);

    % Adjust axes line colors and thickness
    ax = gca;
    ax.XAxis.Color = 'black';
    ax.YAxis.Color = 'black';
    ax.XAxis.LineWidth = 1;
    ax.YAxis.LineWidth = 1;

    % Set tick direction
    ax.XAxis.TickDirection = 'both';
    ax.YAxis.TickDirection = 'both';

    % Add a black border around the plot area
    set(gca, 'Box', 'on', 'LineWidth', 1, 'Color', 'none');

    % Format the data series with pleasant colors
    for j = 1:numel(h)
        set(h(j), 'FaceColor', pleasantColors(j, :), 'EdgeColor', 'none');
    end

    % Add text labels to the bars in scientific notation
    for j = 1:length(FCCs)
        text(FCCs(j) * 1.05, j, sprintf('%.2e', FCCs(j)), 'VerticalAlignment', 'middle', 'FontSize', 12);
    end

    figureName = ['FCC_',ecModel.id,'_',modelFCC.target{i},'.svg'];

    % Save the chart as an SVG file
    saveas(gcf, fullfile(subDirPath,figureName), 'svg');
end