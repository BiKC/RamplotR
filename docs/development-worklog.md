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
- Active branch: `feature/atlas-geometry-state-candidates` from current
  `main` after PR #50. Implements an initial **exploratory** experimental
  geometry comparison, not biological state assignment:
  - extract first-model C-alpha coordinates for exact SIFTS polymer positions
    from the already retrieved PDBe updated mmCIF;
  - retain label asym/seq mapping and ignore ambiguous atoms/positions;
  - require >=30 common observed canonical positions and >=60% common-core
    coverage for every selected structure;
  - rigid-body-invariant C-alpha internal distance-map RMSD on one shared,
    deterministically sampled core of <=300 positions;
  - average-linkage geometric groups with a displayed user-selected distance
    cutoff and representative medoids;
  - pure-R and JS tests, plus real Shiny browser tests with synthetic changed
    internal geometry.
  - Next: scientifically benchmark cluster cutoffs and actual structures,
    inspect construct equivalence and metadata before calling these states.
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

1. Paginated experimental cohort inventory has merged in PR #48.
2. Verify exact PDBe SIFTS residue correspondence for selected entities
   (current active branch), then extend to cohort-level mapping.
   Validate observable residue coverage and accession/isoform consistency
   before treating structures as state evidence.
   Prefer updated PDBe mmCIF/SIFTS residue-level annotations over guessing
   interior mappings from range endpoints.
3. Integrate a source-defined experimental conformational-state clustering
   method and representative structures, with method/provenance recorded.
4. Add fragment-level backbone switch regions and prediction-state clustering.
5. Compare prediction-state coverage against experimentally observed states.
6. Benchmark coarse backbone-state labels before treating them as scientific
   conformation definitions.
7. Update manuscript only after analysis, validation and case studies are
   settled.

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
