# First large-structure scaling results

**Date:** September 28, 2026  
**Workflow:** [Large-structure scaling benchmark, successful run 36430181464](https://github.com/BiKC/RamplotR/actions/runs/36430181464)  
**Dataset:** 6VXX, parsed once using Bio3D and saved to RDS; bundled `original` density reference.

The source has 23,694 atom rows, 2,916 protein residues and 2,844 residues
with classifiable backbone angles. To test scaling, we duplicated the same
source with unique chain names, without changing any atomic coordinates.
All replicated copies produced identical torsions and region labels.

## Recorded measurements

All multipliers ran in separate R processes on **one Ubuntu 24.04 GitHub-hosted
runner**. A process performed three iterations. The first built the reference
profile cache; the other two reused it.

| Copies | Atom rows | Residues | Classified | Cold total (s) | Warm total, both runs (s) | Peak RSS (KiB) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 23,694 | 2,916 | 2,844 | 0.355 | 0.058, 0.054 | 135,628 |
| 3 | 71,082 | 8,748 | 8,532 | 0.400 | 0.090, 0.089 | 141,304 |
| 10 | 236,940 | 29,160 | 28,440 | 0.679 | 0.311, 0.284 | 195,956 |

The classification stage took 0.240, 0.226 and 0.285 seconds on the cold
iterations at 1x, 3x and 10x, respectively. On the two warm iterations,
classification took 0.011–0.012 seconds at 1x, 0.014–0.015 seconds at 3x
and 0.018–0.020 seconds at 10x. The timed backbone stage ranged from
0.042–0.047 seconds (1x warm) to 0.266–0.291 seconds (10x warm).

Every value above is from the attached CSV/log artifact in the workflow run,
not an estimate. Peak RSS is measured by GNU `time` and includes the entire
R process. These measurements exclude the earlier source download and PDB
parsing, and the repeated copies are synthetic. The CI runner is shared-host
infrastructure, so the figures are a **pilot scaling experiment**, not
publication-grade comparisons to the old version or other programs.

## Next measurements for the preprint

- Repeat the same protocol on a dedicated, specified machine, with more
  repetitions and the exact R/Bio3D versions archived.
- Include genuinely independent large structures, not only replicated 6VXX.
- Compare the historical version only after deciding which calculations
  are scientifically comparable across changed classification methods.
- Obtain an independent structural-validation report for the region labels.
