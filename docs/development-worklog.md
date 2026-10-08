# RamplotR development worklog

> Read this file before continuing roadmap implementation after a reconnect or lost chat context.
> Update it whenever a feature is merged, abandoned, or a new blocker is discovered.

## Current baseline

- Current `main` baseline when this work started: `6d1e4e95cdf5df4ae63bdbeccb4a6eeb7af41a1a`.
- Manuscript work remains separate on `arxiv-preprint` / PR #20.
- The roadmap deliberately focuses RamplotR on backbone conformation, predicted-model confidence and comparative structural analysis rather than recreating a full MolProbity/Phenix validation suite.
- PR #43 — Compare workflow progressive-disclosure polish — is merged.
- PR #44 — post-merge mobile Compare regression audit — is merged and green; no CSS changes were required.
- PR #45 — coarse backbone-state transitions — is merged.
- Manuscript PR #20 remains open separately.

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
- PR #43 — integrated Compare workflow polish.
  - comparison source/provenance uses progressive disclosure and auto-collapses after load;
  - Swap remains next to chain-role controls;
  - summary evidence is grouped into Alignment, Backbone, Validation and conditional Prediction confidence.

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

## Active work

- Baseline is current `main` after the canonical-coordinate implementation
  (commit `c2615795ff202e5019d9abc03b130125e77a565b`).
- Canonical mapping (Phase 0A) already exists in `main`:
  `shinyRam/R/canonical.R`, `canonical-mapping.js`, inspector/summary
  integration and tests. Do not reimplement it on the older
  `feature/conformational-atlas-foundation` branch.
- Active branch: `feature/atlas-experimental-inventory`.
- Current slice: Phase 1 inventory discovery for a UniProt accession:
  - experimental polymer-entity discovery from RCSB Search;
  - bounded metadata retrieval, entity/reference coverage and limitations;
  - an Atlas tab with comparison handoff to the existing Compare workflow;
  - unit tests plus scientific CI/UI contract checks.
- This is inventory discovery, not experimental conformational-state clustering.
  It must not report a "state count" until a clustering analysis is implemented.
- Long-term roadmap: `docs/conformational-atlas-roadmap.md`.
- Manuscript remains separate on `arxiv-preprint` / PR #20.
- Before continuing after a reconnect, inspect `main`, recent PRs,
  and this worklog. Do not restore older branches over new work.

## Next implementation steps

1. Finish and verify the experimental inventory branch; merge after the
   scientific regression suite passes.
2. Add exact SIFTS canonical coverage per returned entity/chain, preserving
   unmapped/ambiguous residues and isoform differences.
3. Add pagination and caching so proteins with hundreds of PDB entities are
   analysed as complete, auditable cohorts rather than truncated search hits.
4. Implement state discovery via established archive clusters when available,
   with explicit structural/alignment validation and representative structures.
5. Implement fragment-level switch regions, prediction-state clustering and
   prediction-to-experimental state coverage.
6. Benchmark coarse backbone-state bins before making biological claims.
7. Keep the manuscript branch separate until analysis and validation stabilize.

## Branch hygiene

Before starting or continuing work:

1. read this worklog;
2. inspect open and recently merged PRs;
3. inspect `main` recent commits;
4. search existing branch names for the intended feature;
5. branch from current `main`, unless an intentionally stacked PR is explicitly documented here.

Old feature/fix/chore branches are retained in GitHub but are **not active work** unless this file explicitly says otherwise. This includes, among others:

- `feature/rama8000-validation`
- `feature/conformational-change-explorer`
- `feature/prediction-ensemble` and `feature/prediction-ensemble-reporting`
- `feature/experimental-counterparts` and `feature/experimental-counterpart-search`
- `feature/group-conformation-comparison*`
- `feature/comparison-background-swap` and `fix/compare-background-swap`
- `feature/compare-local-context`, `feature/compare-alignment-quality`,
  `feature/compare-prediction-confidence`, and `feature/compare-workflow-polish`
- `feature/af3-prediction-ensembles`, `feature/local-structure-context`,
  `feature/residue-evidence-inspector`, and `feature/ensemble-evidence-inspector`
- `chore/browser-test-stability-worklog` and `chore/prediction-upload-order`
- `rebuild/compare-prediction-confidence-20261007`

The only long-lived non-main branch intentionally kept for active content is
`arxiv-preprint` / PR #20. The active implementation branch is the one named in **Active work**
above; all older feature branches are historical unless explicitly reopened.
