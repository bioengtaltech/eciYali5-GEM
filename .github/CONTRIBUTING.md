# Contributor Guidelines

First of all, thank you for contributing to **eciYali5-GEM**!  
Anybody is welcome to contribute, but please follow these guidelines.

You can contribute in two main ways:
1. By creating **issues**, and  
2. By sending **pull requests (PRs)** with additions, deletions, or corrections to the model.

---

## Reporting Issues

Please report an issue at  
 https://github.com/bioengtaltech/eciYali5-GEM/issues  
if you notice any of the following:

- Incorrect or missing annotation in the model.  
- Missing feature or field that should be included.  
- Unexpected or inconsistent simulation results.  
- Incomplete or unclear documentation.  
- Any other feedback or bug.

When creating an issue, please make sure that:

- You have tested your code (if applicable) and dependencies.  
- You are using the `main` branch of the repository.  
- You provide any necessary files or links for understanding the issue.  
- You checked that a similar issue does not already exist.  

Feel free to comment on any https://github.com/bioengtaltech/eciYali5-GEM/issues.  
Please also comply with the https://github.com/bioengtaltech/eciYali5-GEM/blob/main/.github/CODE_OF_CONDUCT.md.

If you find this model useful, please consider **starring** the repository — it helps others discover it and lets us know who’s using **eciYali5-GEM**!

---

## Contributing to the Model

If you want to contribute additions or improvements, start by creating an issue and assigning it to yourself to describe what you plan to do.  
This helps avoid duplicated efforts and allows discussion before coding.

You can also check https://github.com/bioengtaltech/eciYali5-GEM/issues for tasks tagged with **help wanted** — these are open for community contribution.

---

### Setting Up the Repository Locally

1. Make sure you have all [requirements](https://github.com/bioengtaltech/eciYali5-GEM#requirements) for running and editing the model.

2. Fork this repository:
   ```
   https://github.com/bioengtaltech/eciYali5-GEM
   ```

3. Clone your fork locally:
   ```bash
   git clone https://github.com/<your-github-username>/eciYali5-GEM.git
   ```

4. Check out the branch you want to base your work on (usually `develop`):
   ```bash
   git checkout develop
   ```

5. Create a new branch for your work:
   ```bash
   git checkout -b feat/my-feature
   ```

6. Make your changes locally.  
   - Model edits should be scripted in MATLAB or Python, **not by directly editing the XML/YAML/TXT files**.  
   - Scripts for curation go in `code/modelCuration/`. You can use `TEMPLATEcuration.m` or adapt a previous `vX_X_X.m` version script.  
   - Place general scripts under `/code`, and tabular data (as `.tsv` files) under `/data`.  
   - Avoid storing binary files (`.mat`, `.xls`, etc.) in the repository.  

7. Before committing:
   - Review your changes with `git diff` or your Git client.  
   - Make sure changes are intentional and clean.  

8. Commit and push your branch:
   ```bash
   git add .
   git commit -m "feat-rxn: add methanol pathway"
   git push origin feat/my-feature
   ```

9. Create a **pull request** (PR) to merge your branch into `develop`:
   - Go to your fork on GitHub.
   - Click **New Pull Request** → choose your branch → base it on `bioengtaltech/eciYali5-GEM:develop`.
   - In the description, summarize your change and link related issues (e.g. *Closes #15*).

---

## Branching Model

| Branch | Purpose |
|:--------|:---------|
| `develop` | The main working branch for ongoing work. |
| `main` | Contains the latest tested and released version. |
| `{chore, doc, feat, fix, refactor, style}/descriptive-name` | For all new work. Examples: `feat/new_pathway`, `fix/reaction_annotations`. |

---

## ✍️ Semantic Commits

Use **concise, descriptive** commit messages following this format:

```
action-object: short description
```

### Common Actions
| Action | Meaning |
|:--------|:---------|
| `feat` | Add a new feature (reaction, metabolite, etc.) |
| `fix` | Correct something wrong in the model |
| `refactor` | Restructure without changing functionality |
| `style` | Format or syntax change only |
| `doc` | Update or add documentation |
| `chore` | Maintenance or dependency update |
| `test` | Add or fix tests |

### Optional Object Tags
`rxn`, `rxn.annot`, `rxn.prop`, `met`, `met.annot`, `met.prop`, `gene`, `gene.annot`, `comp`, `comp.annot`, `data`

### Examples

| Commit type | Example message |
|:-------------|:----------------|
| Add new pathway | `feat-rxn: methanol pathway` |
| Fix duplicated metabolite | `fix-met: duplicated citrate` |
| Update annotations | `fix-gene.annot: updated IDs from UniProt` |
| Add metabolite formulas | `feat-met.prop: carbohydrate formulas` |
| Update toolbox | `chore: update RAVEN version` |

A more detailed explanation can go in the commit description if needed.

---

## Acknowledgments

These contribution guidelines were adapted from  
[SysBioChalmers/yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM),  
with inspiration from  
[opencobra/cobrapy](https://github.com/opencobra/cobrapy) and  
[SysBioChalmers/RAVEN](https://github.com/SysBioChalmers/RAVEN).

---
