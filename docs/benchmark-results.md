# Reproducible benchmark and validation notes

These are measured GitHub Actions results, not projected performance. The original reference distribution was used. Workflow artifacts include exact CSV files and R session details.

## First baseline: PDB 1CRN

Run: https://github.com/BiKC/RamplotR/actions/runs/36423887036
Date: 2026-09-28
Structure: PDB 1CRN (46 residues, 327 atoms as parsed)
Bio3D package installed by runner: 2.4-5

The angle comparison matched 37 phi and 36 psi angles. Median absolute difference was 0 degrees for both. The maximum phi difference was 2.84e-14 degrees and maximum psi difference was 2.84e-13 degrees. Comparison tolerance was 0.5 degrees.

| Iteration | Backbone extraction (s) | Classification (s) | Total (s) |
| --- | ---: | ---: | ---: |
| 1 | 0.140 | 0.221 | 0.361 |
| 2 | 0.048 | 0.171 | 0.219 |
| 3 | 0.049 | 0.189 | 0.238 |

Parsing the downloaded structure took 0.352 seconds. The object occupied 54,376 bytes in R, which is **not** a peak RAM measurement. The workflow classified 44 residues. Timings come from a shared-hosted CI runner and are not stable benchmark comparisons against other software.

## Next dataset

The separate workflow also runs PDB 6VXX, a three-chain SARS-CoV-2 spike structure, for a larger real-data baseline. This is a different virus and protein from the published monkeypox application; do not conflate the datasets.

## Repeat locally

Install Bio3D in a recent R 4.x environment and, from the repository root, run:

```sh
Rscript benchmarks/validate_vs_bio3d.R 1CRN benchmarks/output/local-validation.csv
Rscript benchmarks/run.R 1CRN 3 original benchmarks/output/local-timing.csv
Rscript benchmarks/validate_vs_bio3d.R 6VXX benchmarks/output/viral-validation.csv
Rscript benchmarks/run.R 6VXX 3 original benchmarks/output/viral-timing.csv
```

The benchmark stores CSV timings and the session details needed to identify the R and package versions. Further work is needed for MolProbity region-level comparisons, large-complex scaling and peak RSS measurements.
