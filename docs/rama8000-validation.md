# Rama8000 standard validation

RamplotR evaluates backbone geometry in two deliberately separate ways:

1. **Native RamplotR density regions** use the selected RamplotR reference
   dataset and report Favoured, Allowed, Generously allowed and Not allowed.
2. **Rama8000 standard validation** uses the six current cctbx/Phenix
   Rama8000 residue classes and reports Favored, Allowed and Outlier.

The standard-validation result is independent of the selected RamplotR plotting
background. A native RamplotR `Not allowed` residue and a Rama8000 `Outlier`
are therefore not treated as synonyms.

## Implementation

The Rama8000 implementation mirrors the current cctbx/Phenix code pinned in
`shinyRam/static/rama8000/SOURCE.md`.

The six classes are:

- General
- Gly
- cis-Pro
- trans-Pro
- pre-Pro
- Ile/Val

Class assignment follows cctbx priority. Proline is classified as cis-Pro when
the preceding peptide omega angle is strictly between -90 and +90 degrees;
otherwise it is trans-Pro. Glycine is handled separately. A non-Gly/non-Pro
residue immediately before a bonded proline is pre-Pro; remaining Ile/Val
residues use their shared class; all others are General.

The six vendored 180 x 180 Rama8000 score tables are sampled every two degrees
at odd-numbered phi/psi coordinates from -179 to +179 degrees. RamplotR uses
the same periodic bilinear interpolation as cctbx.

Current score thresholds are reproduced directly:

| Rama8000 class | Favored | Allowed | Outlier |
| --- | ---: | ---: | ---: |
| General | score >= 0.0200 | 0.0005 <= score < 0.0200 | score < 0.0005 |
| cis-Pro | score >= 0.0200 | 0.0020 <= score < 0.0200 | score < 0.0020 |
| Gly, trans-Pro, pre-Pro, Ile/Val | score >= 0.0200 | 0.0010 <= score < 0.0200 | score < 0.0010 |

## Direct comparison with official wwPDB reports

The existing independent validation corpus was rerun with the Rama8000
implementation on 6 October 2026. Coordinate files and official wwPDB
validation XML are independently downloaded and matched by model, chain,
residue number, insertion code and residue identity.

| Structure | Method | Comparable residues | Matching Rama8000/wwPDB categories | Agreement |
| --- | --- | ---: | ---: | ---: |
| 1CRN | X-ray | 44 | 44 | 100% |
| 1UBQ | X-ray | 74 | 74 | 100% |
| 6VXX | cryo-EM | 2,844 | 2,844 | 100% |
| 2DQ4 | X-ray | 682 | 682 | 100% |
| 1D3Z | NMR, model 1 | 74 | 74 | 100% |
| **Total** |  | **3,718** | **3,718** | **100%** |

The same 3,718 residues also reproduce the independently reported phi/psi
angles to within 0.055 degrees. The category comparison is now direct:
no mapping of the four native RamplotR density regions is involved.

The outlier-containing 2DQ4 sample provides the clearest regression case.
All seven official wwPDB Ramachandran outliers are independently classified
as Rama8000 Outlier by RamplotR. Their native RamplotR density regions remain
visible separately (six Generously allowed, one Allowed), demonstrating why
the two outputs must not be conflated.

The exact standard-validation snapshot is pinned in
`validation/rama8000-baseline-2026-10-06.csv`. The workflow fails if the
same pinned source files no longer reproduce 100% category agreement.

## Reproduce

From the repository root:

~~~bash
Rscript -e 'install.packages(c("bio3d", "xml2", "digest"))'
Rscript benchmarks/compare_wwpdb.R 2DQ4 original benchmarks/output/wwpdb/2DQ4
Rscript tests/wwpdb-baseline.R 2DQ4 benchmarks/output/wwpdb/2DQ4
~~~

Outputs include:

- `rama8000_comparison.csv`: residue-wise standard category, class, score and
  official wwPDB category;
- `rama8000_contingency.csv`: direct Favored/Allowed/Outlier contingency;
- `rama8000-summary.csv`: coverage and exact category agreement;
- the existing coordinate/validation source hashes and session metadata.

## Interpretation

Rama8000 is the standard validation result intended for statements about
Favored/Allowed/Outlier status. The native RamplotR density regions remain
useful for exploring alternative bundled density distributions and for
RamplotR-specific visualization, but `Not allowed` is not renamed or
interpreted as a Rama8000 outlier.
