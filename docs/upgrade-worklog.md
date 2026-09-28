# RamplotR upgrade worklog

This log records the main changes and their verification status. Work happens in separately scoped commits. The published historical baseline is preserved as `v0.1.0-legacy` (commit `aa3eba2180d505f4dc01c9c2bf977cef00b3252a`).

## Completed in PR #2

1. Residue-aware classification independent of the plotting background, with an explicit legacy option.
2. Deterministic density thresholds and reusable reference grids.
3. Atom-based backbone torsions that respect peptide connectivity, chain boundaries and insertion codes.

## Phase 2: CI, inputs, reproducibility

- [x] Keep scientific changes in their own commits, merged without squashing.
- [x] Add a two-platform R regression workflow for existing scientific test scripts.
- [ ] Verify the workflow succeeds on both supported CI operating systems.
- [x] Add targeted tests for classification with synthetic reference grids and edge cases.
- [x] Detect local PDB versus mmCIF input based on its original uploaded filename.
- [x] Add clear upload errors and an explicit PDB-versus-upload source choice.
- [x] Update installation instructions and list direct R dependencies; a reproducible renv lockfile is pending a real R environment.
- [ ] Record independent scientific comparison and benchmark results before the preprint.

The scientific classifications are tied to the bundled reference densities. Do not claim MolProbity-equivalent outlier percentages without a separate benchmark. New methods and reference datasets will be versioned.

## Phase 3: scientific comparisons and timing

- [x] Provide a repeatable analysis/timing script recording the structure, reference dataset, software environment and per-stage elapsed time.
- [x] Provide an independent torsion-angle comparison with Bio3D for PDB accessions without insertion codes.
- [ ] Run the comparison on multiple public proteins and investigate mismatches.
- [ ] Run scaling benchmarks on large PDB/mmCIF inputs and record benchmark output.
- [ ] Compare region labels against independently generated structural-validation reports; numerical thresholds differ across implementations.
