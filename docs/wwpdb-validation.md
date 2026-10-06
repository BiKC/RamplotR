# Independent wwPDB Ramachandran validation

RamplotR compares its calculations against official wwPDB validation XML using
a pinned five-structure corpus. The workflow keeps three questions separate:

1. are phi/psi angles calculated correctly?
2. does RamplotR's Rama8000 implementation reproduce current standard
   Favored/Allowed/Outlier categories?
3. how do the optional native RamplotR density regions relate to the same
   residues?

Only the first two are scientific validation gates. Native density regions are
reported independently and are never remapped into an artificial outlier class.

## Public reproducible corpus

The committed `validation/manifest.csv` identifies:

- 1CRN and 1UBQ — X-ray crystallography;
- 6VXX — a large cryo-EM spike complex;
- 2DQ4 — X-ray crystallography with seven official Ramachandran outliers;
- 1D3Z — solution NMR, model 1.

The GitHub Actions workflow downloads each RCSB mmCIF and the corresponding
official wwPDB validation XML. Every run preserves source files, SHA256 hashes,
download URLs, git revision and R package versions.

Residues are joined by model, chain, residue number, insertion code and residue
identity. Alternate-conformation handling is conservative; ambiguous matches
reduce coverage instead of being silently paired.

## Outputs

For each accession the workflow produces:

- `residue_comparison.csv`: independent angle comparison plus native
  RamplotR regions;
- `rama8000_comparison.csv`: direct Rama8000 versus wwPDB categories;
- `rama8000_contingency.csv`: direct Favored/Allowed/Outlier contingency;
- `summary.csv`: angle coverage and maximum circular differences;
- `rama8000-summary.csv`: standard-validation coverage and category agreement;
- `provenance.csv` and `r-session.txt`: exact source and software provenance.

## Gates

Angle validation requires at least 90% coverage of finite RamplotR phi/psi
pairs and no circular difference above 1.5°. The pinned corpus currently has
100% coverage and a maximum observed difference of 0.055°.

Rama8000 validation requires **100% category agreement** with the official
wwPDB report for the pinned corpus. On 6 October 2026 this is 3,718/3,718
comparable residues.

Source hashes and the quantitative snapshots are additionally checked by
`tests/wwpdb-baseline.R`. A changed official source must be investigated and
reviewed before any baseline is updated.

## Reproduce

~~~bash
Rscript -e 'install.packages(c("bio3d", "xml2", "digest"))'
Rscript tests/wwpdb.R
Rscript benchmarks/compare_wwpdb.R 1CRN original benchmarks/output/wwpdb/1CRN
Rscript tests/wwpdb-baseline.R 1CRN benchmarks/output/wwpdb/1CRN
~~~

For exact replay, pass the pinned compressed XML and mmCIF as the fourth and
fifth arguments to `benchmarks/compare_wwpdb.R`.

## Scope

The corpus is a functional validation set, not a population sample of the PDB.
Only model 1 is used in the independent benchmark. This validates the current
angle calculation and Rama8000 implementation; it does not make RamplotR a
replacement for broader MolProbity/wwPDB validation of clashes, rotamers,
covalent geometry or experimental data quality.

See [the current results](validation-results.md) and
[the Rama8000 implementation details](rama8000-validation.md).
