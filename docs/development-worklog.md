# RamplotR development worklog

> Read this file before continuing roadmap implementation after a reconnect or lost chat context.
> Update it whenever a feature is merged, abandoned, or a new blocker is discovered.

## Current baseline

- Current `main` baseline: `860ee1bb95f5af336db42528be7ef81aa2d3db53`.
- Manuscript work remains separate on `arxiv-preprint` / PR #20.
- The roadmap deliberately focuses RamplotR on backbone conformation, predicted-model confidence and comparative structural analysis rather than recreating a full MolProbity/Phenix validation suite.
- As of 7 October 2026 there are **no open implementation PRs**; only the manuscript PR remains open.

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
- PR #33 — local ligand/hetero-residue structural context in the main inspector.
- PR #34 — AlphaFold 3 sample ensemble analysis.
- PR #36 — comparison Ramachandran density background and quick structure-role swap.
- PR #37 — prediction-ensemble disagreement linked into the residue inspector.
- PR #38 — compare nearby ligand/hetero context across aligned residues.
- PR #39 — pairwise sequence identity and chain-coverage reporting with interpretation warnings.
- PR #40 — compare prediction confidence across aligned residues.
  - Explicit prediction provenance for uploaded comparison structures.
  - pLDDT on both sides, signed delta-pLDDT, confidence-change filtering and plot/table/inspector integration.

Do **not** reimplement these features on old feature branches; several merged branches still exist.

## Browser-test blocker: resolved

The shared prediction-upload timeout was fixed centrally in PRs #41 and #42.

Current stable browser-test ordering:

1. choose prediction provenance first;
2. wait for the reactive upload controls to settle;
3. attach the file;
4. verify the native `File` object exists;
5. allow Shiny to finish the upload before submitting.

PR #38 was rebuilt against this central fix and passed its scientific and browser workflows before merge.
PR #40 was rebuilt cleanly on current `main`; scientific, structure-validation and full live-browser workflows passed on the final feature code. The independent wwPDB workflow had also passed on the prior identical scientific feature implementation; the final rerun was still waiting in runner setup when #40 was merged.

## Next implementation steps

1. **Audit/polish the integrated Compare workflow.**
   - It now combines source/provenance, chain alignment quality, conformational-change track, matched density background, 2D/3D linked selection, Rama8000, pLDDT differences, ligand/hetero context, filters/table and group comparison.
   - Prefer progressive disclosure and clearer information hierarchy over adding more permanent panels.
2. Verify the post-merge `main` browser experience at normal laptop width and mobile width after any Compare polish.
3. Only start a new scientific feature after checking open PRs, recent main commits and branch names.
4. Synchronize the manuscript only after the implementation/UI settles.
   - Describe current Rama8000 standard validation, not old/native outlier interpretations.
   - Position RamplotR around residue-centred backbone comparison, predicted-model confidence and linked structural context.

## Branch hygiene

Before starting or continuing work:

1. read this worklog;
2. inspect open and recently merged PRs;
3. inspect `main` recent commits;
4. search existing branch names for the intended feature;
5. branch from current `main`, unless an intentionally stacked PR is explicitly documented here.

Old merged branches such as `feature/rama8000-validation`,
`feature/conformational-change-explorer` and `feature/prediction-ensemble`
are historical implementation branches and must not be treated as active work.

A temporary branch named `rebuild/compare-prediction-confidence-20261007`
was used to reconstruct PR #40 cleanly on top of the then-current `main`.
It is not active development work.
