# RamplotR development worklog

> Read this file before continuing roadmap implementation after a reconnect or lost chat context.
> Update it whenever a feature is merged, abandoned, or a new blocker is discovered.

## Current baseline

- Main branch baseline when this worklog was created: `9b4f77f200888f779208a9c249900abe865d4af0`.
- Manuscript work remains separate on `arxiv-preprint` / PR #20.
- The roadmap intentionally focuses RamplotR on backbone conformation, predicted-model confidence and comparative structural analysis rather than recreating a full MolProbity/Phenix validation suite.

## Roadmap features already merged

- PR #21 — six-class Rama8000 standard validation.
  - General, Gly, cis-Pro, trans-Pro, pre-Pro and Ile/Val.
  - Direct cctbx/Phenix-style Favored/Allowed/Outlier evaluation.
  - Pinned five-structure validation reproduces 3,718/3,718 wwPDB categories.
- PR #22 — Conformational Change Explorer.
  - Wrapped delta-phi/delta-psi, combined local backbone displacement, clickable alignment track and linked 2D/3D focus.
- PR #23 + #25 — multi-file prediction ensembles and standalone reporting.
- PR #28 — discovery of experimental PDB counterparts for predicted models.
- PR #30 — residue evidence inspector.
- PR #32 — structure-group conformational comparison.
- PR #33 — local ligand/hetero-residue structural context.
- PR #34 — AlphaFold 3 sample ensemble analysis.
- PR #36 — comparison Ramachandran density background and quick structure-role swap.
- PR #37 — prediction-ensemble disagreement linked into the residue inspector.

Do **not** reimplement these features on old feature branches; several old branches still exist after merge.

## Currently open feature PRs

- PR #38 — compare local ligand/hetero context across aligned residues.
  - Scientific tests pass.
  - Browser preview currently fails in a shared file-upload wait in `tests/prediction-browser.cjs`; not currently evidence of a PR-specific scientific bug.
- PR #39 — show pairwise alignment identity and coverage.
  - Scientific, structure-validation and wwPDB workflows pass.
  - Same shared browser-test upload wait fails.
- PR #40 — compare prediction confidence across aligned residues.
  - Scientific, structure-validation and wwPDB workflows pass.
  - Same shared browser-test upload wait fails.

## Shared blocker discovered 2026-10-07

All three open PRs time out at `tests/prediction-browser.cjs` line 95 while waiting for Shiny's internal file-upload progress/input representation after Puppeteer `uploadFile()`.

The file is already present in the native browser input, so the test should not depend solely on Shiny's private `$inputValues["structfile:shiny.file"]` key or progress-bar state. Fix this once on main/test infrastructure, then rerun/rebase feature PRs rather than patching each feature independently.

## Next implementation steps

1. Stabilize the prediction browser file-upload wait centrally.
2. Rerun PRs #38, #39 and #40; fix only feature-specific failures that remain.
3. Merge completed PRs when all workflows pass.
4. Review the integrated Compare workflow for information density and duplicated controls after #38-#40 are combined.
5. Continue roadmap only after checking merged/open branches to avoid duplicating work.
6. Keep the manuscript synchronized only after implementation settles; do not describe old/native validation outputs as if they were the current standard Rama8000 result.

## Branch hygiene

Before starting a new feature:
1. inspect open and recently merged PRs;
2. search branch names for the intended feature;
3. read this worklog;
4. branch from current `main`, not from an old merged feature branch unless intentionally stacking work.
