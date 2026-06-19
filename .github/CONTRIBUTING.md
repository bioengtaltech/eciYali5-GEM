# Contributor Guidelines

Thank you for contributing to **eciYali5-GEM**. Contributions are welcome
through issues and pull requests (PRs) with additions, deletions, corrections,
or documentation improvements.

## Reporting Issues

Please report issues at https://github.com/bioengtaltech/eciYali5-GEM/issues
if you notice any of the following:

- Incorrect or missing model annotation.
- Missing model feature or field.
- Unexpected or inconsistent simulation results.
- Incomplete or unclear documentation.
- Any other feedback or bug.

When creating an issue, please make sure that:

- You have tested your code and dependencies when applicable.
- You are using the `main` branch when reporting behavior in the latest release.
- You provide files, links, commands, or scripts needed to understand the issue.
- You checked that a similar issue does not already exist.

If you find this model useful, please consider **starring** the repository — it helps others discover it and lets us know who’s using **eciYali5-GEM**!
Please comply with the
[code of conduct](https://github.com/bioengtaltech/eciYali5-GEM/blob/main/.github/CODE_OF_CONDUCT.md)
when participating in issues, PRs, and reviews.

## Contributing To The Model

If you want to contribute additions or improvements, start by creating an issue
and assigning it to yourself to describe what you plan to do. This helps avoid
duplicated effort and allows discussion before coding.

You can also check https://github.com/bioengtaltech/eciYali5-GEM/issues for
tasks tagged with **help wanted**.

## Setting Up The Repository Locally

1. Make sure you have the
   [required software](https://github.com/bioengtaltech/eciYali5-GEM#required-software)
   for running and editing the model.

2. Fork this repository:

   ```bash
   https://github.com/bioengtaltech/eciYali5-GEM
   ```

3. Clone your fork locally:

   ```bash
   git clone https://github.com/<your-github-username>/eciYali5-GEM.git
   ```

4. Check out `develop`, the base branch for ongoing work:

   ```bash
   git checkout develop
   ```

5. Create a branch for your work:

   ```bash
   git checkout -b feat/my-feature
   ```

6. Make your changes locally.

   - Model edits must be scripted in MATLAB or Python, not made by directly
     editing generated XML/YML/TXT model files.
   - Curations to released models go in `code/modelCuration/`. Copy
     `TEMPLATEcuration.m` to a versioned script such as `v1_0_1.m`, then run it
     with `runCurationScript`.
   - General reusable MATLAB helpers go in `code/utilities/`.
   - Publication-specific scripts and outputs currently live under
     `code/debiaggi2025/`.
   - Tabular data should be stored as `.tsv` or `.csv` files under `data/` or,
     for one-off curation support files, under `data/modelCuration/`.
   - Do not store generated model `.mat` or `.xlsx` files on `develop`.
     Source model files under `model/**/sourceModels/`, such as
     `model/iYali5-GEM/sourceModels/iYali4_corr.mat`, may remain when they are
     intentionally required by the generation pipeline.

7. Review your changes before committing:

   ```bash
   git diff
   ```

8. Commit and push your branch:

   ```bash
   git add .
   git commit -m "feat(rxn): add methanol pathway"
   git push origin feat/my-feature
   ```

9. Create a PR to merge your branch into `develop`:

   - Go to your fork on GitHub.
   - Click **New Pull Request** and choose
     `bioengtaltech/eciYali5-GEM:develop` as the base branch.
   - Summarize your change and link related issues, for example `Closes #15`.

## Curation Workflow

This repository includes a lightweight Yeast-GEM-style curation scaffold:

```matlab
addpath(genpath('code'))
model = runCurationScript('v1_0_1', 'ModelName', 'eciYali5-GEM', 'Save', true);
```

The curation runner loads a model, executes the selected script, runs
`validateEciYaliModel`, checks the `develop` binary-file policy, and can save
updated XML/YML/TXT model files through `saveEciYaliModel`.

Use `getEarlierModelVersion` when a curation needs to start from a tagged model
without switching the current branch.

## Branching Model

| Branch | Purpose |
|:-------|:--------|
| `develop` | Base branch for ongoing work and PRs. |
| `main` | Latest tested and released version. |
| `{chore,docs,feat,fix,refactor,style,test}/descriptive-name` | Branches for new work, for example `feat/new-pathway` or `fix/reaction-annotations`. |

## Conventional Commits

Use Conventional Commits 1.0.0 for commit messages:

```text
<type>[optional scope]: <description>
```

Preferred types are `feat`, `fix`, `docs`, `test`, `refactor`, `chore`,
`build`, `ci`, `perf`, `style`, and `revert`.

Examples:

| Change | Example commit message |
|:-------|:-----------------------|
| Add a pathway | `feat(rxn): add methanol pathway` |
| Fix a duplicated metabolite | `fix(met): remove duplicated citrate` |
| Update gene annotations | `fix(gene.annot): update IDs from UniProt` |
| Add metabolite formulas | `feat(met.prop): add carbohydrate formulas` |
| Update contributor docs | `docs(contributing): clarify curation workflow` |
| Update toolbox metadata | `chore(deps): update RAVEN version` |

Use a commit body when more detail, references, or issue links are useful.

## Releases

Documentation and workflow-helper changes do not require a new GEM release by
themselves. A new model release is needed when model content, release artifacts,
or `version.txt` intentionally change and maintainers decide the updates are
ready for `main`.

## Acknowledgments

These contribution guidelines were adapted from
[SysBioChalmers/yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM), with
inspiration from [opencobra/cobrapy](https://github.com/opencobra/cobrapy) and
[SysBioChalmers/RAVEN](https://github.com/SysBioChalmers/RAVEN).
