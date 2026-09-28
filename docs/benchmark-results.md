# Reproducible benchmark and validation notes

These are measured GitHub Actions results, not projected performance. The original reference distribution was used. Workflow artifacts include exact CSV files and R session details.

## Corrected, coverage-checked Bio3D validation (2026-09-28)

The first comparison script joined residue identifiers as exact strings.
Bio3D left-pads some residue identifiers with spaces, so early results
undercounted comparable positions. The corrected script trims only surrounding
whitespace, preserves insertion-code distinctions, rejects duplicate keys and
requires at least 98% coverage of eligible residue IDs and finite RamplotR
angles. These are the corrected results from
[GitHub Actions run 36428648286](https://github.com/BiKC/RamplotR/actions/runs/36428648286);
[final cross-platform checks and structural validation](https://github.com/BiKC/RamplotR/actions/runs/36428693471)
also passed.

| Structure | RamplotR residues | Matched residue IDs | Matched phi / own finite | Matched psi / own finite | Maximum phi difference (°) | Maximum psi difference (°) | Differences over 0.5° |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1CRN | 46 | 46 | 45 / 45 | 45 / 45 | 2.84e-14 | 2.84e-13 | 0 |
| 6VXX | 2916 | 2916 | 2880 / 2880 | 2880 / 2880 | 1.42e-13 | 1.14e-11 | 0 |

Bio3D reports 2,979 residue rows for 6VXX, while RamplotR extracts
2,916 supported amino-acid residues. The coverage percentages above
therefore refer to the RamplotR-eligible residue IDs and finite angles,
not to every row in Bio3D's output. This comparison verifies torsion
angles. It does not validate reference-density choices or region labels
against independent structural-validation software.

### Timings from the corrected validation run

| Structure | Parsing, measured once (s) | Iteration | Backbone (s) | Classification (s) | Processing total (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1CRN | 0.260 | 1, cold profile | 0.053 | 0.168 | 0.221 |
| 1CRN | 0.260 | 2, cached | 0.002 | 0.008 | 0.010 |
| 1CRN | 0.260 | 3, cached | 0.002 | 0.007 | 0.009 |
| 6VXX | 0.550 | 1, cold profile | 0.067 | 0.162 | 0.229 |
| 6VXX | 0.550 | 2, cached | 0.019 | 0.008 | 0.027 |
| 6VXX | 0.550 | 3, cached | 0.018 | 0.008 | 0.026 |

Parsing takes place once before the three iterations; its timing is
repeated for context in the exported CSV. The processing total excludes
parsing. Reference RDS files are warmed before the iterations, but the
reference-classification profile is first built in iteration 1. These
measurements come from a shared GitHub-hosted runner and are not a
controlled cross-version performance comparison. Object size is not peak
resident memory.

## Historical first baseline: PDB 1CRN (partially matched identifiers)

The figures in this historical section were produced before the whitespace
correction and must not be reported as full-coverage validation.

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
