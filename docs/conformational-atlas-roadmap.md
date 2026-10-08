# Conformational Atlas roadmap

## Product direction

RamplotR should remain a focused backbone-conformation workbench rather than
becoming another general validation, refinement or molecular-visualisation
suite.

The long-term target is a third research workflow alongside single-structure
inspection and direct comparison:

> Given a protein, reconstruct the experimentally observed and predicted
> conformational states, identify the local backbone regions that distinguish
> those states, connect the states to available structural context, and show
> whether a prediction ensemble captures or misses the experimentally observed
> alternatives.

The central object remains the residue or short backbone segment. Global
structure information is added only when it helps explain those local
conformational states.

## Progress

- Phase 0A canonical residue mapping is implemented on `main`. It provides
  conservative PDBe SIFTS/UniProt mapping and direct AlphaFold DB coordinates.
- Phase 1 initial experimental inventory was merged in PR #47. It fetches
  experimental polymer entities by UniProt accession, reports
  method/resolution and reference sequence coverage, and opens a chosen
  candidate in Compare.
- Experimental inventory pagination is implemented on
  `feature/atlas-pagination-cohort` and is pending CI/merge. It preserves
  page offsets, deduplicates polymer entities, records failed metadata IDs,
  detects changes to the archive hit count, and loads the next page on demand.
- **Not implemented yet:** exact per-entity canonical SIFTS residue coverage,
  experimental structural-state clustering, fragment switch regions and
  prediction-vs-experiment coverage.
- The manuscript remains separate until those scientific methods are
  implemented and validated.

## Scientific principles

1. Keep evidence types separate. Rama8000 validation, native RamplotR density,
   prediction confidence, experimental metadata and conformational variation
   must never be collapsed into one opaque quality score.
2. Never interpret prediction-model counts as thermodynamic populations.
3. Never call an experimentally unobserved predicted state wrong solely because
   it has not been deposited in the PDB.
4. Do not infer causation from structural metadata. Ligand, mutation or partner
   enrichment is association unless supported by a study design that justifies
   a stronger statement.
5. Preserve PDB author numbering and insertion codes while adding canonical
   UniProt coordinates. Do not replace one numbering system with the other.
6. Do not fabricate residue mappings. SIFTS range mappings are expanded only
   where a one-to-one author-number/UniProt-number relation is provable.
7. Coarse Alpha-R/Beta/PPII/Alpha-L states remain human-readable navigation
   bins until separately benchmarked. They are not secondary-structure or
   validation assignments.
8. Archive-scale or ensemble clustering must expose coverage, unresolved
   residues and alignment uncertainty.

## Phase 0: foundations

### 0A. Canonical residue coordinates

Goal: every structure used in an atlas can refer to the same protein coordinate
system.

- Add a canonical mapping data model alongside PDB chain/residue/insertion IDs.
- Add PDBe SIFTS UniProt mapping retrieval for PDB accessions.
- Prefer the current PDBe v2 mapping endpoint and keep a documented compatibility
  fallback while necessary.
- Normalize mapping segments with accession, entity, author chain,
  struct-asym ID, UniProt range, PDB/author range, sequence identity and
  coverage.
- Expand segment mappings to exact residue mappings only when that expansion is
  unambiguous.
- Add direct UniProt coordinates for AlphaFold DB models, whose coordinate
  residue numbers follow the requested UniProt model sequence.
- Show canonical coordinates in the residue inspector when available.
- Preserve mapping provenance and coverage in exports.
- Later add exact residue-level SIFTS parsing from PDBe-enriched mmCIF/SIFTS
  data for segments that cannot safely be expanded.

Acceptance criteria:

- PDB numbering and insertion codes remain unchanged.
- Canonical mappings never silently shift residues.
- Ambiguous/nonlinear ranges remain explicitly unresolved.
- Mapping helpers are pure and unit tested.
- Browser/Shinylive retrieval does not require a new server dependency.

### 0B. Benchmark the coarse backbone-state layer

Goal: keep the current state-change UI scientifically useful without turning a
simple navigation scheme into an unsupported biological classification.

- Build a benchmark set spanning common and unusual backbone regions.
- Quantify stability of the current nearest-centre labels around state borders.
- Compare the coarse labels against a fragment-level structural alphabet.
- Keep the current labels as overview/navigation bins unless the benchmark
  supports a stronger interpretation.
- Record the versioned centres and cutoff in provenance.

## Phase 1: experimental Conformational Atlas

Goal: start from a UniProt accession or mapped loaded structure and collect the
experimental structural evidence automatically.

- Add an Atlas input based on UniProt accession.
- Enumerate relevant experimental PDB polymer entities.
- Retrieve canonical residue coverage and basic entry/entity metadata.
- Retain method, resolution, construct/mutation information, ligand/cofactor
  context, assembly/partner context and source identifiers where available.
- Cache fetched metadata within the session.
- Make archive/network use explicit to the user.
- Do not add a permanent Atlas tab until this workflow can return a useful
  result rather than a placeholder.

Initial result:

- number of experimental entries/entities;
- canonical sequence coverage;
- methods/resolution overview;
- mapped structures suitable for state analysis;
- structures excluded because of insufficient mapping/alignment.

## Phase 2: experimental state discovery

Goal: identify distinct experimentally observed conformational states without
reimplementing archive infrastructure unnecessarily.

- Investigate and reuse PDBe/PDBe-KB conformational/superposition clusters when
  available for the target.
- Retain the cluster method and source provenance.
- Provide a local fallback only when archive clusters are unavailable.
- Select representative structures for each state.
- Calculate state coverage on the canonical UniProt coordinate system.
- Show whether differences are mainly local remodeling, global/domain movement,
  or both.

Outputs:

- experimental state count;
- structures supporting each state;
- representative structure;
- canonical coverage;
- global separation metric;
- local switch regions.

## Phase 3: local conformational fingerprints and switch regions

Goal: move beyond isolated single-residue angle changes.

- Keep the coarse Alpha-R/Beta/PPII/Alpha-L/Other state as an overview layer.
- Add a fragment-level torsion representation, initially targeting a
  five-residue structural alphabet or an equivalent explicitly benchmarked
  torsion fingerprint.
- Detect contiguous regions whose local backbone fingerprints differ between
  states.
- Report state consensus and missing-data coverage per region.
- Link each switch region back to the existing Ramachandran, sequence and 3D
  inspector.

Example output:

    residues 101-108
    experimental state A: local fingerprint X in 18/20 structures
    experimental state B: local fingerprint Y in 14/16 structures

## Phase 4: prediction-ensemble state discovery

Goal: answer how many distinct structural states a prediction workflow sampled.

- Superpose compatible prediction models using a documented structure-aware
  method.
- Build a model-to-model structural distance representation.
- Cluster prediction models into structural states.
- Keep residue-level phi/psi spread and confidence as separate evidence.
- Choose representative prediction models per state.
- Link cluster-defining switch regions to the existing residue inspector.
- Support AF2/ColabFold, AF3 sample sets, ESMFold and AFsample-style ensembles.
- Report model counts as sample counts, never populations.

## Phase 5: prediction versus experiment state coverage

Goal: determine whether prediction ensembles reproduce experimentally observed
states.

- Match predicted states to experimental states using transparent structural
  criteria.
- Report captured experimental states.
- Report experimentally observed states not sampled by the prediction ensemble.
- Report predicted states that are currently unobserved experimentally without
  calling them incorrect.
- Add residue/switch-region-level agreement so missed states can be localized.

Primary research output:

    Experimental state 1: captured
    Experimental state 2: captured
    Experimental state 3: not sampled
    Prediction state P3: no close experimental counterpart

## Phase 6: structural-context association

Goal: explain what available evidence is associated with state differences.

Candidate context:

- ligand/cofactor presence;
- protein/nucleic-acid interaction partner;
- interface membership;
- engineered mutation;
- construct boundaries;
- experimental method and resolution;
- assembly/oligomeric context;
- deposited condition metadata when sufficiently standardized.

Analysis:

- descriptive contingency tables first;
- effect sizes and confidence intervals where assumptions are satisfied;
- Fisher exact tests for suitable categorical comparisons;
- circular/statistical analyses only when sample size and independence permit.

Wording must remain association-based unless causality is externally
established.

## Phase 7: global versus local movement

Goal: cover the main blind spot of a phi/psi-centred analysis.

- Add a global/domain-rearrangement signal, preferably based on
  transformation-independent inter-residue distances.
- Identify hinge-like regions when possible.
- Separate:
  - global/domain rearrangement with little local remodeling;
  - local backbone remodeling;
  - mixed changes.

This layer complements local torsions rather than replacing them.

## Phase 8: AF3 and complex-aware states

Goal: extend the same conformation-first workflow to interaction context.

- Compare the same protein across monomer/complex predictions.
- Add interface membership and partner-contact context to switch regions.
- Retain intra/inter-chain confidence evidence separately from local
  stereochemistry.
- Compare predicted interfaces with experimental interface annotations when
  canonical mapping permits it.

## Phase 9: reproducibility, automation and publication

- Atlas CSV/JSON/HTML exports with complete provenance.
- Offline/scriptable Atlas analysis for reproducible cohorts.
- Persistent cache format with source timestamps/checksums.
- Benchmarks for archive-scale targets and large prediction ensembles.
- Browser tests for the end-to-end Atlas workflow.
- Curated biological case studies where the inferred state differences are
  already well established experimentally.
- Update the manuscript only after the workflow and benchmark claims are
  stable.

## Deliberate non-goals

Do not turn RamplotR into:

- a full MolProbity/Phenix replacement;
- a refinement or model-building program;
- a generic MD trajectory analyser;
- a generic multiple-structure alignment package;
- a generic PDB/UniProt annotation browser;
- a variant-pathogenicity predictor;
- a general-purpose molecular viewer.

External tools and archive services should be reused when they already solve
those problems well.

## Implementation order

1. Canonical SIFTS/UniProt mapping.
2. Coarse-state benchmark.
3. Experimental Atlas retrieval and metadata normalization.
4. Experimental state discovery.
5. Fragment-level switch-region analysis.
6. Prediction-state clustering.
7. Prediction-to-experiment state coverage.
8. State/context association.
9. Global/domain-motion context.
10. AF3/interface-aware extensions.
11. Full reporting, batch mode, biological case studies and manuscript update.

Each phase should be merged independently with tests and an updated development
worklog. Old feature branches must not be reused as active work.
