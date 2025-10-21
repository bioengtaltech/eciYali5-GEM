%% STAGE 1: expansion from a starting metabolic model to an ecModel structure
% Preamble
binaryFiles = true;

% STEP 1 Set modelAdapter
adapterLocation = fullfile(pwd,'eciYali5GEMadapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

% STEP 2 Load conventional iYali
model = loadConventionalGEM();

% To avoid unecerray DLKcat reruns
DLKcat = false;

% STEP 3-4 Prepare ecModel. NB: I generated a custom uniprot.tsv file where I
% got the KEGG crossref and exchanged it for the gene_oln
% pseudoRxns.tsv file is also manually made
[ecModel, noUniprot] = makeEcModel(model,false);

% STEP 5 Store model in YAML format
%saveEcModel(ecModel,'eciYali_stage1.yml');

%% STAGE 2: integration of kcat into the ecModel structure

% STEP 6 Fuzzy matching with BRENDA
% Requires EC numbers, which are here first taken from the starting model,
% with the missing ones taken from Uniprot & KEGG databases.
ecModel         = getECfromDatabase(ecModel);

% Do the actual fuzzy matching with BRENDA.
kcatList_fuzzy  = fuzzyKcatMatching(ecModel);

% STEP 7 DLKcat prediction through machine learning
% Requires metabolite SMILES, which are gathered from PubChem.
[ecModel, noSmiles] = findMetSmiles(ecModel);

if DLKcat == true
    % DLKcat runs in Python. An input file is written, which is then used by
    % DLKcat, while the output file is read back into MATLAB.
    writeDLKcatInput(ecModel,[],[],[],[],true);

    % runDLKcat will run the DLKcat algorithm via a Docker image. If the
    % DLKcat.tsv file already has kcat values, these will all be overwritten.
    runDLKcat();
end

kcatList_DLKcat = readDLKcatOutput(ecModel);

% STEP 8 Combine kcat from BRENDA and DLKcat
kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy, 2); % Changed DLKcat priority

% STEP 9 Take kcatList and populate edModel.ec.kcat
ecModel  = selectKcatValue(ecModel, kcatList_merged);
ecModel  = selectKcatValue(ecModel, kcatList_DLKcat,'max','ifHigher'); % another way of prioritizing DLKcat

% STEP 10 Apply custom kcat values
[ecModel, rxnUpdated, notMatch] = applyCustomKcats(ecModel);

% STEP 11 Get kcat values across isozymes
ecModel = getKcatAcrossIsozymes(ecModel);

% STEP 12 Get standard kcat
[ecModel, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(ecModel);

% STEP 13 Apply kcat constraints from ecModel.ec.kcat to ecModel.S
ecModel = applyKcatConstraints(ecModel);

% STEP 14 Set upper bound of protein pool
% The protein pool exchange is constrained by the total protein content
% (Ptot), multiplied by the f-factor (ratio of enzymes/proteins) and the
% sigma-factor (how saturated enzymes are on average: how close to their
% Vmax to they function based on e.g. metabolite concentrations). In
% modelAdapter Ptot, f- and sigma-factors can all be specified (as rough
% estimates, 0.5 for each of the three parameters is reasonable).
Ptot  = params.Ptot;
f     = params.f;
sigma = params.sigma;

ecModel = setProtPoolSize(ecModel,Ptot,f,sigma);

%% STAGE 3: model tuning

% STEP 15 Test maximum growth rate
% Test whether the model is able to reach maximum growth if glucose uptake
% is unlimited. First set glycerol uptake unconstraint.
ecModel = setParam(ecModel,'lb','1714',0);
ecModel = setParam(ecModel,'lb','1808',-1000);

% And set growth maximization as the objective function.
ecModel = setParam(ecModel,'obj','xBIOMASS',1);
% Run FBA.
sol = solveLP(ecModel,1);
bioRxnIdx = getIndexes(ecModel,params.bioRxn,'rxns');
fprintf('Growth rate: %f /hour.\n', sol.x(bioRxnIdx))

% STEP 16 Sensitivity tuning
ecModel = setProtPoolSize(ecModel);
[ecModel, tunedKcats] = sensitivityTuning(ecModel);

% Inspect the tunedKcats structure in table format.
struct2table(tunedKcats)

%% Save model files
% Fetch this script's directory
scriptDir = fileparts(mfilename('fullpath'));

% Go two folders up (from /code/debiaggi2025/ to /eciYali5-GEM/)
baseDir = fileparts(fileparts(scriptDir));

% Create path for iYali5-GEM
ecGEMpath = fullfile(baseDir, 'model');

% Save files
if binaryFiles == false
    exportForGit(model,'eciYali5-GEM',ecGEMpath,{'xml','yml','txt'},false,false);
else
    exportForGit(model,'eciYali5-GEM',ecGEMpath,{'xml','yml','txt','xlsx','mat'},false,false);
end