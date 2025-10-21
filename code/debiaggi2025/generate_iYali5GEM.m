% generate_iYali5GEM.m
% Update script for the iYali4_corr genome-scale metabolic model (GEM)
% Targets: prepare model for GECKO compatibility / RAVEN-based workflows

clear
clc

% 
binaryFiles = false;

% Fetch this script's directory
scriptDir = fileparts(mfilename('fullpath'));

% Go two folders up (from /code/debiaggi2025/ to /eciYali5-GEM/)
baseDir = fileparts(fileparts(scriptDir));

% Build the path to the model directory
modelsDirPath = fullfile(baseDir, 'model', 'iYali5-GEM', 'sourceModels');

% COBRA format fixes — load and normalise imported model
modelFile = fullfile(modelsDirPath, 'iYali4_corr.mat');
load(modelFile)  % local .mat copy to avoid import issues

% Load iYali4 (Kerkhoven, 2016) model
modelOG = readCbModel(fullfile(modelsDirPath,'iYali4.xml')); % original COBRA xml import

% Correct data type of model.subSystems (cell array of strings required)
var = modelOG.subSystems;
for i = 1:length(model.subSystems)
    var{i} = model.subSystems(i);
end
model.subSystems = var;

% Append compartment code to metabolite IDs to avoid ambiguous IDs across compartments
for i = 1:numel(model.mets)
    model.mets{i} = [model.mets{i}, model.Comps{i}];
end

% Create and correct compartment metadata
model.comps = modelOG.comps;
model.compNames = modelOG.compNames;

% Map selected compartment indices to standard short codes (assumption: indices align as intended)
model.comps([1, 2, 3, 5, 6, 9, 11, 12, 13]) = {'cy', 'en', 'ex', 'em', 'go', 'mi', 'nu', 'pe', 'va'};

%% RAVEN format fixes
% Convert COBRA model to RAVEN-compatible structure
model = ravenCobraWrapper(model);

% Overwrite specific compartments with single-letter codes used in RAVEN (c: cytosol, e: extracellular, m: mitochondrion)
model.comps([1, 3, 9]) = {'c', 'e', 'm'};

%% Add biomass pseudometabolites and pseudoreaction
% Add carbohydrate pseudoreaction
model = addMetabolites(model, 'm1836', 'carbohydrate', 'c', 'Pseudometabolite for biomass equation. Necessary for sumBiomMass function of RAVEN based GEMs');
model = addReaction(model, 'xCARBOHYDRATE', {'m294','m401','m1123','m1324','m1836'}, [-0.00686, -0.868358, -0.943397, -0.235849, 1], 'carbohydrate pseudoreaction');

% Add DNA pseudoreaction
model = addMetabolites(model, 'm1837', 'DNA', 'c', 'Pseudometabolite for biomass equation. Necessary for sumBiomMass function of RAVEN based GEMs');
model = addReaction(model, 'xDNA', {'m89','m459','m465','m505','m1837'}, [-0.01007, -0.010436, -0.009226, -0.010353, 1], 'DNA pseudoreaction');

% Add RNA pseudoreaction
model = addMetabolites(model, 'm1838', 'RNA', 'c', 'Pseudometabolite for biomass equation. Necessary for sumBiomMass function of RAVEN based GEMs');
model = addReaction(model, 'xRNA', {'m86','m93','m95','m149','m1838'}, [-0.055401, -0.050871, -0.05721, -0.059363, 1], 'RNA pseudoreaction');

% Add ion pseudoreaction
model = addMetabolites(model, 'm1839', 'ion', 'c', 'Pseudometabolite for biomass equation. Necessary for sumBiomMass function of RAVEN based GEMs');
model = addReaction(model, 'xION', {'m964','m1839'}, [-0.02, 1], 'ion pseudoreaction');

% Modify protein pseudoreaction
rxns = {'xAMINOACID'};
equations.mets = {'m50','m74','m114','m130','m267','m272','m310','m319','m441','m443',...
    'm615','m743','m765','m770','m772','m775','m793','m859','m992','m1008','m1726'};
equations.stoichCoeffs = {[-0.243904, -0.044232, -0.567939, -0.243866, -0.186531, -0.509943, -0.125563,...
    -0.186498, -0.284699, -0.003632, -0.053881, -0.090566, -0.206372, -0.210543,...
    -0.002154, -0.19455, -0.275005, -0.082571, -0.041282, -0.172773, 1]};

model = changeRxns(model, rxns, equations, 1);

% Change the name to protein
% get xAMINOACID index
idx = find(strcmp(model.rxns,rxns));

model.rxns(idx) = {'xPROTEIN'};

model.rxnNames(idx) = {'protein pseudoreaction'};

% Modify lipid pseudoreaction
rxns = {'xLIPID'};
equations.mets = {'m1000','m1579','m1581','m1640','m1641','m1648','m1651','m1705','m1700','m1701','m1631','m1727'};
equations.stoichCoeffs = {[-0.021176,-0.0001,-0.000305,-0.01251,-0.003152,-0.0036,-0.000842,-0.003399,-0.010101,-0.01424,-0.001967,1]};

model = changeRxns(model, rxns, equations, 1);

% Change the name to lipid
% get xLIPID index
idx = find(strcmp(model.rxns,rxns));

model.rxnNames(idx) = {'lipid pseudoreaction'};

% Remove old lipid pseudoreaction
model = removeReactions(model,'2108',false,false,false);

% Modify xBIOMASS
rxns = {'xBIOMASS'};
equations.mets = {'m32', 'm141', 'm1836', 'm1837', 'm1838', 'm1839', 'm1726', 'm1727',...
    'm10', 'm35', 'm143', 'm1401'};
equations.stoichCoeffs = {[-23.09, -23.09, -1, -1, -1, -1, -1, -1, 23.09, 23.09, 23.09, 1]};

model = changeRxns(model, rxns, equations, 1);

% Change the name to biomass
% get xBIOMASS index
idx = find(strcmp(model.rxns,rxns));

model.rxnNames(idx) = {'biomass pseudoreaction'};

%% Reaction restorations and additions

% Malic enzyme (NADPH) was improperly removed. Restore the original
% reaction and fix the name of malic enzyme (NADH). For more details, see
% (https://doi.org/10.1007/s10529-013-1302-7)
idx = find(strcmp(model.rxns,'718'));
model.rxnNames(idx) = {'malic enzyme (NAD)'};
model = addReaction(model, '719', {'m538','m178','m06','m44','m176'}, [-1, -1, 1, 1, 1], 'malic enzyme (NADP)');
model = changeGeneAssoc(model,'719','YALI0E18634g');

% Glycerol-3-phosphate dehydrogenase (FAD) - UniID Q6CEQ0 - is annotated as cytosolic as well
model = addReaction(model, 'G3PD', {'m528','m562','m456','m1048'}, [-1, -1, 1, 1], 'glycerol-3-phosphate dehydrogenase (FAD)');
model = changeGeneAssoc(model,'G3PD','YALI0B13970g');

% Adds lipid exchange reaction so it can be maximized by FSEOF. subSystem
% must be fixed manually
[model, addedRxns] = addExchangeRxns(model,'out','m1640');
model.subSystems(end) = model.subSystems(end-1);

% Remove rxns
% There is no evidence that Yarrowia encodes a acyl
% dihydroxyacetonephosphate reductase. There was no GPR for this rxn, and
% this enzyme cannot be found in Yarrowia's Uniprot.
model = removeReactions(model,'iYL0336',false,false,false);

% Methylglyoxal synthase is not found in Yarrowia (Uniprot search). This
% reaction is a secondary reaction of TPI but in eucaryotes it seems to be
% present only in the animal kingdom.
model = removeReactions(model,'1936',false,false,false);

% Old biomass equations
model = removeReactions(model,{'2110','2133','4041','biomass_C','newBiom'},false,false,false);

% Reversibility changes
% The reversibility of alcohol dehydrogenases must be dealt with
% separate reactions
model = setParam(model,'rev',{'163','165'},0);

% Fixing the reversibility of isocitrate dehydrogenase (NAD+)
model = setParam(model,'rev','658',0);

%% Fix rxn bounds
% DAG acyltransferase rxn (lb) is not reversible
model = setParam(model,'lb',{'336u'},0);
model = setParam(model,'lb',{'2948g','2884g'},0);

% DAG acyltransferase rxn (mm) was blocked for no apparent reason
model = setParam(model,'ub',{'336u_1','337u_2'},[1000,1000]);

% Block glucose uptake because we are working with glycerol
model = setParam(model,'eq','1714',0);

% Allow glycerol uptake
model = setParam(model,'lb','1808', -1000);
model = setParam(model,'ub','1808', 1000);

% Glycerol dehydrogenase. It is not native to Yarrowia.
model = setParam(model,'eq','487', 0);

% Cytosolic aconitase. GPR that points to it reveals a mitochondrial
% enzyme in Uniprot
model = setParam(model,'eq','303', 0);

% Cytosolic aconitate hydratase. GPR that points to it reveals a mitochondrial
% enzyme in Uniprot
model = setParam(model,'eq','2305', 0);

% Aldehyde dehydrogenases. The cytosolic ones are probably unused.
model = setParam(model,'eq',{'173','2116'}, 0);

% Isocitratre dehydrogenase is only mitochondrial in Yali
% (https://doi.org/10.1007/s12010-013-0373-1). The cytosolic and
% peroxissomal reactions are YeastGEM artifacts.
model = setParam(model,'eq',{'659','661'}, 0);

% See reversibility tab
model = setParam(model,'lb',{'163','165','658'}, 0);

% Were blocked for no apparent reason.
model = setParam(model,'ub',{'495','2512g','275u','277u','3789g','8u'}, [1000,1000,1000,1000,1000,1000]);

% Reaction names fixes
model.rxnNames(find(strcmp(model.rxns,'27'))) = {'homoaconitase'};
model.rxnNames(find(strcmp(model.rxns,'117'))) = {'2-methylcitrate dehydratase'};
model.rxnNames(find(strcmp(model.rxns,'iYL0459'))) = {'phosphoribosylglycinamide formyltransferase 1'};

% GPR fixes
model = changeGeneAssoc(model,'958','YALI0C24101g'); % had another GPR associated to an enzyme with a different function
model = changeGeneAssoc(model,'336u','YALI0E32769g or YALI0D07986g'); % GPR for DGA2, which is missing in the original model
model = changeGeneAssoc(model,'336u_1','YALI0E32769g or YALI0D07986g'); % GPR for DGA2, which is missing in the original model
model = changeGeneAssoc(model,'yli0053','YALI0C14806g');
model = changeGeneAssoc(model,'2117','YALI0C05258g or YALI0E20977g'); % Both proteins catalyse this reactions (https://doi.org/10.1111/1751-7915.13745)
model = changeGeneAssoc(model,'2119','YALI0C05258g or YALI0E20977g'); % Both proteins catalyse this reactions (https://doi.org/10.1111/1751-7915.13745)
model = changeGeneAssoc(model,'923','YALI0F09559g'); % Removed enzyme from the complex

% Ethanol dehydrogenases
% The reverse rxn (acetaldehyde to ethanol is not as ubiquitous as the
% forward. GPR kept only for the protein with the greatest homology to
% S. cerevisae ADH1, shown experimentally to be the only ADH of yeast
% to metabolize acetaldehyde
% (https://doi.org/10.1111/j.1567-1364.2011.00760.x)
model = changeGeneAssoc(model,'2115','YALI0A16379g');

% Mitochondrial alcohol dehydrogenases
% Some genes (YALI0D01738g, YALI0C06171g, YALI0E12463g, YALI0A15147g, YALI0F25003g, YALI0F08129g, YALI0B08404g, YALI0E19921g) were removed due to low homology according to Uniprot
model = changeGeneAssoc(model,'165','YALI0F29623g');

% Aldehyde dehydrogenases
rxns = {'172';...
    '173';...
    '174';...
    '175';...
    '201';...
    '2116'};

grRules = {'YALI0D07942g or YALI0F04444g';...
    'YALI0C03025g';...
    'YALI0E00264g';...
    'YALI0E00264g';...
    'YALI0D07942g or YALI0F04444g';...
    'YALI0D07942g or YALI0F04444g';...
    };

model = changeGrRules(model,rxns,grRules);

% Phosphoribosylglycinamide formyltransferase 1 - iYL0459
% A0A1D8ND08 must be found in the uniprot.tsv file and the gene
% YALI1D03865g added to it.
model = changeGeneAssoc(model,'iYL0459','YALI1D03865g');

% Set growth as objective funcion
model = setParam(model,'obj',{'xBIOMASS'},1);

% Calculate GAMnonPol
GAMnonPol = calculateGAM(model);

% Fix biomass composition
[X,~,C,~,~,~,~] = sumBioMassYali4(model, false);
delta = X - 1;  % difference to balance
fC = (C - delta) / C;
model = rescalePseudoReaction(model, 'carbohydrate', fC);

% Recalculate and set GAM
[~,~,model] = calculateGAM(model, GAMnonPol,true);

% Test changes
sol = solveLP(model, 1);
model = setParam(model, 'obj', 'xBIOMASS', 1);
sol2 = solveLP(model, 1);

% Name model
model.modelID = 'iYali5';
model.description = 'iYali5.xml';
model.modelName = 'v0.1';

% Create path for iYali5-GEM
GEMpath = fullfile(baseDir, 'model', 'iYali5-GEM');

% Save files
if binaryFiles == false
    exportForGit(model,'iYali5',GEMpath,{'xml','yml','txt'},false,false);
else
    exportForGit(model,'iYali5',GEMpath,{'xml','yml','txt','xlsx','mat'},false,false);
end
%% Helper function to add metabolites
function model = addMetabolites(model, metID, metName, compartment, metNotes)
metsToAdd.mets = {metID};
metsToAdd.metNames = {metName};
metsToAdd.compartments = {compartment};
metsToAdd.metNotes = {metNotes};

model = addMets(model, metsToAdd);
end

% Helper function to add reactions
function model = addReaction(model, rxnID, mets, stoichCoeffs, rxnName)
rxnsToAdd.rxns = rxnID;
rxnsToAdd.mets = mets;
rxnsToAdd.stoichCoeffs = stoichCoeffs;
rxnsToAdd.rxnNames = {rxnName};
rxnsToAdd.lb = 0;
rxnsToAdd.ub = 1000;
rxnsToAdd.subSystems = {''};

model = addRxns(model, rxnsToAdd);
end