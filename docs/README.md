# RamplotR documentation

Choose a guide by the research task you want to perform. The [main README](../README.md) introduces the application, screenshots and installation.

## Working with structures

- [Interactive inspection and publication-ready figures](inspection-user-guide.md): linked Ramachandran plot, residue table, sequence navigator, 3D viewer, comparison and exports.
- [AlphaFold, ColabFold and ESMFold confidence](prediction-confidence.md): model imports, pLDDT and linked PAE inspection.
- [Structural verification, cryo-EM overlays and ensembles](structural-verification.md): measured backbone/side-chain geometry, independent official reports, maps and offline batch processing.

## Scientific validation and performance

- [Independent wwPDB validation method](wwpdb-validation.md) and [five-structure results](validation-results.md).
- [Benchmark results](benchmark-results.md), [scaling study](scaling-study.md) and [scaling measurements](scaling-results.md).

## Development records

[Implementation and verification history](development/README.md) is kept separately from the user guides. Its historical notes describe what changed during development, not different current application modes.

- [Atlas cluster sensitivity](atlas-cluster-robustness.md): deterministic position-deletion checks, shared-PDB warnings, CSV evidence and interpretation limits.

- [Atlas experimental context](atlas-experimental-context.md): recorded method, resolution, chemistry and same-PDB caveats without invented ligand or functional states.
