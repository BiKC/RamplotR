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

- PR #46: canonical SIFTS/UniProt coordinates merged.
- PR #47: first experimental Conformational Atlas inventory merged.
  - UniProt-to-experimental-polymer-entity search;
  - RCSB method, resolution and entity/reference coverage metadata;
  - bounded, capped retrieval with explicit partial-result notices;
  - comparison handoff and desktop/mobile regression coverage.
- PR #48: experimental Atlas cohort pagination merged into `main`
  (`d70ba566fa2330e0d238413d803e880109a7c20d`).
  - Paged RCSB search, deterministic sort/offset tracking, deduplicated entity
    IDs, archive-total consistency checks and metadata failure identifiers.
  - Browser, scientific (Ubuntu/Windows), independent wwPDB and structure
    validation workflows all passed. The first browser attempt had a transient
    prediction-upload timeout; rerun passed the full suite including Atlas.
- PR #49 merged into `main` (`daa014eba9f92eab8cd79cb32466cc5970e2af1c`):
  exact PDBe updated-mmCIF SIFTS mapping for selected experimental entities,
  with insertion-safe joins, conflict reporting and on-demand verification.
  Scientific Ubuntu/Windows, structure, independent wwPDB and live-browser CI
  all passed.
- PR #50 merged into `main` (`67935cba2fdcf954fbc4afc133cdeefa54e98c90`).
  Aggregates verified SIFTS mappings by UniProt position, retains conflicts
  in exact residue exports, excludes ambiguous coordinates from observed
  position-support counts and provides two CSV exports.
  Scientific Ubuntu/Windows, structure/wwPDB, and the full live-browser Atlas
  regression all passed after a browser test selector fix.
- PR #51 merged into `main`
  (`869b1af85e31b0b770e50b25e8c71ef37f288971`):
  - extracted first-model C-alpha coordinates in PDBe updated mmCIF for
    exact SIFTS polymer positions (no author-number interpolation);
  - required >=30 common observed UniProt positions and >=60% core coverage;
  - rigid-body-invariant C-alpha distance-map RMSD on a consistent sampled core;
  - exploratory average-linkage geometric groups with an explicit adjustable
    threshold and representative entity for each group;
  - two-structure distance plot and >=3-structure dendrogram;
  - Ubuntu/Windows scientific tests, independent wwPDB, structure benchmark
    and live Shiny browser workflow all passed after fixing two plot regressions.
- PR #52 merged to `main` (`8c721a7067b6062720cbf417d1e07474848f5805`):
  - model-1 N/CA/C extraction alongside exact SIFTS-mapped C-alpha;
  - φ/ψ calculation guarded by native polymer adjacency, full backbone
    atoms and 1.0–1.9 Å C–N peptide continuity;
  - residue-wise circular Δφ/Δψ at shared UniProt positions between
    the first two exploratory geometric group representatives;
  - adjustable angular-change threshold, candidate contiguous regions,
    canonical residue plot/table and CSV export;
  - explicit missing-data, ambiguity and source-provenance rules.
- Scientific Ubuntu/Windows tests, structure benchmark, independent wwPDB
  validation and full Shiny browser suite all passed on the final PR #52
  commit. CI now runs the Atlas browser regression before the older prediction
  upload test; the latter also passed on the final run.
- PR #53 merged into `main` (`3fe89c94b2e93acb080a7e0ab4d8e3a6352e1585`).
  Real PDBe updated-mmCIF benchmark using E. coli adenylate kinase:
  4AKE (open apo), 1AKE (closed inhibitor-bound), UniProt P69441.
  - 214 exact canonical residues for each structure and chain, using live
    SIFTS mapping and first-model N/CA/C atomic coordinates.
  - 4AKE-A vs 1AKE-A global C-alpha distance-map RMSD **6.5077 Å**;
    within-4AKE A/B **0.4022 Å**, within-1AKE A/B **0.2513 Å**;
    self-control exactly **0 Å**. All compared on 214 shared positions.
  - Full φ/ψ pairs at 212 residues. >=30° local changes: open/closed 37;
    open same-crystal chain controls 16; closed controls 10; self 0.
  - Global distance-map change per canonical position mean:
    LID 8.600 Å, NMP 7.591 Å, CORE 4.540 Å. These are descriptive
    geometry measurements, not local φ/ψ change or calibrated state scores.
  - Live PDBe retrieval and benchmark GitHub CI green; Ubuntu/Windows
    scientific regression CI green.
  - Case study and exact SHA256 source provenance in
    `docs/atlas-adenylate-kinase-benchmark.md`.
  - This is **one known case**, not threshold validation. Chain A/B within
    a crystal are not independent biological replicates. Differences above
    30° already occur in same-crystal chain controls.
- PR #54 merged into `main` (`37e6f23506814a69e415b8c9a445ce104a47e20c`):
  - Curated, provenance-checked **live experimental** benchmark across
    three proteins: ADK (4AKE/1AKE), maltose-binding protein (1OMP/1ANF)
    and ribose-binding protein (1URP/2DRI), using exact PDBe SIFTS mapping.
  - Matched experimental controls: ADK and RBP within-crystal chain
    comparisons, MBP separately crystallized apo 1OMP/1JW4, and three
    numerical identity controls; synthetic mismatched UniProt/construct
    examples must be rejected.
  - Shared aligned UniProt C-alpha coordinates and complete peptide torsions,
    with per-pair denominators, separate global/local measurements and
    residue-name equality checks. The seven structures had 0 observed
    residue-name mismatches at the compared mapped positions.
  - Open/closed vs control C-alpha dRMSD: ADK 6.5077 Å vs 0.4022/0.2513 Å;
    MBP 2.8060 vs 0.4119 Å; RBP 2.7443 vs 0.0773 Å.
  - Local φ/ψ ≥30° candidate counts: ADK 37/212 vs 16/212 and 10/212;
    MBP 23/368 vs 13/368; RBP 10/269 vs 0/269.
  - Reproducible real PDBe source SHA256, atomic data, residue CSVs and
    report archived by GitHub Actions. All three CI workflows passed:
    scientific Ubuntu/Windows, structure validation and live multi-protein
    case-study retrieval/analysis.
  - Detailed results and limitations in
    `docs/atlas-multiprotein-benchmark.md`. **No universal thresholds
    or statistical sensitivity/specificity estimates** have been established.
    Two of the three proteins are related periplasmic binding proteins.
- PR #55 merged into `main`
  (`71de4a4f3300d39254777edf48231ab3101afc38`).
  - User-visible **Guide** tab with seven task-oriented sections: single
    structure/residue inspection, predicted confidence/ensembles, paired
    Compare, grouped comparisons, experimental Atlas, exports and scientific
    interpretation. Entry links from input, Compare and Atlas; tab-navigation
    action links for the workflows.
  - Group Comparison is organized as reference/chain match, two condition
    upload sets, then analysis. Fixed labels/field widths at responsive
    desktop/laptop/mobile sizes. Hide export actions before analysis exists.
  - Also corrected the paired comparison-source button clipping at narrow
    laptop content widths. Initial "New here?" input link disappears once
    a structure is loaded, preserving the compact analysis toolbar.
  - Removed misleading stale Atlas text claiming geometric grouping and
    local candidate detection are not yet implemented.
  - Real headless-browser regression checks no clipping/overlap at 1440,
    900 and 390 px and verifies guide navigation. Screenshots in browser
    artifact from run 37779032286.
  - Full Interface browser preview, scientific Ubuntu/Windows, and structure
    validation workflows passed on final head. Independent wwPDB checks
    1CRN, 1UBQ, 1D3Z and 6VXX passed; 2DQ4 remained in GitHub runner
    R setup at merge time, without a reported code/test failure.
  - No scientific metric, residue classification, or cutoff was changed.
- No active feature PR after #55. Next: link local Atlas backbone candidates
  to selected 3D/Ramachandran residues, verify deployed browser PDBe CORS,
  and continue independent-fold controls for state classification claims.
- IMPORTANT: this repository merge does not update the static Shinylive
  deployment automatically. Regenerate the export and publish its built files
  before expecting the Guide tab at bikc.be/RamplotR.
- Live PDBe API/CORS and updated mmCIF schema still need a deployment smoke
  test. CI mocks the network responses.
- Do not equate verified PDB-entity counts with independent observations or
  distinct conformations. No experimental state clustering exists yet.
- Browser tests mock PDBe responses; confirm live PDBe CORS and updated-mmCIF
  endpoint on deployment before advertising archive-wide state analysis.
- Manual live RCSB API smoke test remains desirable; CI uses mocked RCSB
  search/metadata responses to make network-independent UI regressions.
- Roadmap: `docs/conformational-atlas-roadmap.md`.
- Browser CI for PR #47 passed. Windows scientific tests passed; Ubuntu
  scientific runner for the last CSS/test-only update had not completed setup
  at merge time. An earlier Ubuntu scientific run on the functional changes
  passed.
- Manuscript PR #20 remains separate; it does not yet describe an Atlas state
  discovery algorithm.


## October 8 feature branch

- Active branch: `feature/atlas-residue-inspection-manual-ci`, based on main `0aceaed8`.
- Atlas now offers selecting any two distinct geometry representatives, displays exact SIFTS author residue identifiers and paired angles in a selectable table, and loads a chosen second representative in Compare when a primary structure exists.
- The Compare handoff does not load the first representative into the primary viewer, does not automatically select the aligned canonical position, and requires checking comparability. Full coordinated 3D focus remains unfinished.
- CI policy revised: seven validation/benchmark workflows run for a newly opened non-draft PR or on `ready_for_review`, not on pushes or subsequent PR commits (`synchronize` excluded). Draft PRs skip jobs until marked ready. Manual `workflow_dispatch` remains as fallback. The README screenshot curation workflow is intentionally manual-only because it writes generated artifacts back to the repository. Use PR #56 readiness to trigger final CI and review results before merge.

## Next implementation steps

1. Verify **deployed-browser** PDBe and RCSB CORS in RamplotR's actual
   Shinylive hosting environment; Node/GitHub server retrieval is not enough.
2. Let users select any pair of experimental Atlas geometry groups or their
   representatives; show exact canonical/local residue identifiers and
   coordinated 3D/Ramachandran inspection for candidate switch regions.
3. Expand biological benchmarks to **independent folds** and independent
   experimental same-state replicates, ideally including constructs and
   intentionally mismatched isoforms. Report clear denominators and
   background, not a universal threshold selected from three proteins.
4. Implement cohort-wide mutation, construct and isoform equivalence checks
   before classifying archive structures as comparable biological states.
5. Compare prediction ensembles with validated experimental geometric
   alternatives only when reference-state definitions are defensible.
6. Keep the manuscript separate until the validation supports its claims.

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
`arxiv-preprint` / PR #20. Start any new work from current `main`, update **Active work** with the actual
branch, and do not reuse historical feature branches.
