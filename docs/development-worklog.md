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
- Active branch: `feature/atlas-adk-biological-benchmark`.
  - Real experimental 4AKE open and 1AKE closed adenylate kinase case study,
    UniProt P69441, via the **same browser SIFTS/mapping parser** in Node.
  - Reproducible PDBe updated-mmCIF retrieval with SHA256 provenance and
    model-1 N/CA/C backbone coordinates only.
  - Use explicit observed UniProt position overlap; no author-number
    interpolation. Compare open/closed plus within-4AKE and within-1AKE
    chain copies and a zero-distance self control.
  - Report global distance-map RMSD separately from local wrapped φ/ψ
    changes, above-threshold regions and domain-level descriptive context.
  - Live-archive GitHub Actions workflow uploads reports, data and hashes.
  - Work is **benchmarking, not cutoff calibration**. Intra-crystal copies are
    not independent observations, and one pair cannot establish a general
    validated state classifier.
- After benchmark CI passes, review the actual metrics, limitations and
  whether the live updated-mmCIF schema matches the parser. No scientific
  claims until the real-data job is green.
- Earlier roadmap after PR #52:
  1. real biological case study and numeric calibration with E. coli
     adenylate kinase open 4AKE vs closed 1AKE (literature-supported pair);
  2. construct/isoform compatibility and full-cohort comparison;
  3. arbitrary pair of experimental group medoids and local switch-region
     inspection in linked 3D;
  4. prediction-to-experimental state coverage;
  5. a separately benchmarked fragment-level structural alphabet.
- Current switch regions remain exploratory, based on two representatives
  and a user-selectable navigation threshold, not a statistical or functional
  state classification.
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

## Next implementation steps

1. Validate live PDBe/RCSB endpoints and CORS in the deployed browser;
   browser CI currently mocks external responses.
2. Benchmark experimental geometry and local backbone differences with curated
   real structure pairs, starting with E. coli adenylate kinase 4AKE (open)
   versus 1AKE (closed). Include same-state replicates and construct controls
   before calibrating any distance or angular thresholds.
3. Extend local φ/ψ comparison to user-selected pairs of geometric groups,
   then all supporting structures with uncertainty/coverage reporting.
4. Make switch-region selection navigate to both corresponding 3D structures,
   Ramachandran plots and precise canonical/local residue identifiers.
5. Add prediction state clustering and experimental state coverage only after
   the experimental reference-state method is validated.
6. Keep the manuscript separate until biological benchmark claims are sound.

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
