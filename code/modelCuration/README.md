# Model Curation Pipeline

This folder contains scripts that consolidate model curations between
released versions, following the workflow used by Yeast-GEM and adapted for
eciYali5-GEM.

## Workflow

1. Copy `TEMPLATEcuration.m` to a versioned script such as `v1_0_1.m`.
2. Add curation code to the versioned script. Keep changes scripted; do not
   edit generated `model/*.xml`, `model/*.yml`, or `model/*.txt` files by hand.
3. Store one-off curation data as TSV files under `data/modelCuration/` when
   the data belongs to the model history.
4. Run the script through `runCurationScript` from MATLAB:

```matlab
addpath(genpath('code'))
model = runCurationScript('v1_0_1', 'ModelName', 'eciYali5-GEM', 'Save', true);
```

By default, `saveEciYaliModel` exports XML, YML and TXT files. Use
`IncludeBinary`, `true` only for release workflows that intentionally need
MAT/XLSX files on `main`.

## Curation Targets

Curations should declare which model they target:

- `iYali5`: base GEM curation. Save `model/iYali5-GEM/iYali5.*`, then
  regenerate `eciYali5-GEM` from the updated base model.
- `eciYali5-GEM`: ecGEM-only curation, such as enzyme constraints, kcats,
  protein-pool settings or ec-model annotations.
- `both`: coordinated work where the base GEM and ecGEM outputs are both
  intentionally updated.

Repository Git tags are snapshots of the whole project. Model-level versions
and derivation are tracked in `model/versions.tsv`, so users can tell which
`iYali5` version an `eciYali5-GEM` release was generated from.

Earlier models can be loaded from a tag without switching branches:

```matlab
iModel = getEarlierModelVersion('v1.0.1', 'iYali5');
ecModel = getEarlierModelVersion('v1.0.1', 'eciYali5-GEM');
```

## Helper Functions

- `loadEciYaliModel` loads `eciYali5-GEM` or `iYali5` from repository model files.
- `saveEciYaliModel` exports model files with repository defaults.
- `getEarlierModelVersion` loads a model from a tag or branch using a temporary
  git worktree without switching the current branch.
- `readModelVersions` reads `model/versions.tsv` and can filter by model name.
- `validateEciYaliModel` runs lightweight structural checks and a solver check
  when `solveLP` is available.
- `assertDevelopModelPolicy` fails on `develop` if generated model binaries are
  present outside `model/**/sourceModels`.
