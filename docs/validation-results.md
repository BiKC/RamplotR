# Independent wwPDB validation: current results

The current validation separates **angle correctness**, **Rama8000 standard
validation** and **native RamplotR density regions**.

The five pinned structures are 1CRN, 1UBQ, 6VXX, 2DQ4 and model 1 of 1D3Z.
Exact coordinate and official wwPDB validation-report SHA256 hashes are stored
in [the source manifest](../validation/baseline-source-hashes-2026-09-28.csv).

## Backbone angles

| Structure | Method | Matched finite phi/psi pairs | Maximum circular difference |
| --- | --- | ---: | ---: |
| 1CRN | X-ray | 44/44 | 0.052° |
| 1UBQ | X-ray | 74/74 | 0.054° |
| 6VXX | cryo-EM | 2,844/2,844 | 0.055° |
| 2DQ4 | X-ray | 682/682 | 0.055° |
| 1D3Z | NMR, model 1 | 74/74 | 0.055° |

All **3,718/3,718** eligible pairs reproduce the independently reported wwPDB
angles to within 0.055°. The official XML prints angles to one decimal place.

## Rama8000 categories

RamplotR's six-class Rama8000 implementation is compared directly with the
official wwPDB Favored/Allowed/Outlier labels. No native RamplotR category
mapping is used.

| Structure | Comparable residues | Exact category matches | Agreement |
| --- | ---: | ---: | ---: |
| 1CRN | 44 | 44 | 100% |
| 1UBQ | 74 | 74 | 100% |
| 6VXX | 2,844 | 2,844 | 100% |
| 2DQ4 | 682 | 682 | 100% |
| 1D3Z | 74 | 74 | 100% |
| **Total** | **3,718** | **3,718** | **100%** |

The exact snapshot is pinned in
[rama8000-baseline-2026-10-06.csv](../validation/rama8000-baseline-2026-10-06.csv).
CI requires the same pinned sources to continue reproducing every category.

For 2DQ4, all seven official wwPDB Ramachandran outliers are also classified
as **Rama8000 Outlier** by RamplotR. Their native RamplotR density regions are
shown separately: six are Generously allowed and one is Allowed. The exact
residue identities and both current outputs are pinned in
[2dq4-wwpdb-outlier-examples.csv](../validation/2dq4-wwpdb-outlier-examples.csv).

## Interpretation

Use **Rama8000** when referring to standard Favored/Allowed/Outlier validation.
Use **RamplotR density regions** when discussing the selected exploratory
density reference. The native label `Not allowed` is not a synonym for
Rama8000/wwPDB `Outlier`.

See [Rama8000 standard validation](rama8000-validation.md) for the implementation
and [the reproducible wwPDB protocol](wwpdb-validation.md) for source handling.
