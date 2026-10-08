# Atlas experimental benchmark across three proteins

This benchmark tests whether RamplotR's global distance-map and local
backbone-angle methods behave sensibly across several known experimental
conformational contrasts. **It is not a trained or validated functional-state
classifier.**

## Curated case set

| Protein | Open structure | Closed structure | Available experimental control |
| --- | --- | --- | --- |
| *E. coli* adenylate kinase (P69441) | 4AKE, apo | 1AKE, Ap5A-inhibited | Chain copies within both crystals |
| *E. coli* maltose-binding protein (P0AEX9) | 1OMP, apo | 1ANF, maltose-bound | 1OMP vs independent apo 1JW4 crystal |
| *E. coli* ribose-binding protein (P02925) | 1URP, apo | 2DRI, ribose-bound | 1URP same-crystal chain copies |

One exact self-comparison is added for every protein. It must have zero
geometry difference and zero local-angle changes.

Maltose-binding and ribose-binding protein are distinct proteins but belong
to related periplasmic binding-protein structural families. This sample is
not a representative collection of independent folds.

## Method

A curated, version-controlled JSON manifest fixes the entries, UniProt
accessions, experimental context and control design.

1. Download *updated* mmCIF from PDBe and record each original file's SHA256.
2. Join SIFTS label chain/sequence coordinates to exact UniProt positions.
3. Extract observed model-1 N, Cα and C coordinates.
4. Require at least 100 common, unambiguous canonical Cα positions and 60%
   common-core coverage of each structure.
5. Cross-check deposited monomer identities at common canonical positions.
   Report modifications and substitutions separately, and reject cases with
   more than 5% observed residue-name disagreement.
6. Calculate intrachain Cα distance-map RMSD (Å), matched circular Δφ/Δψ,
   and the number and fraction of comparable positions exceeding the existing
   exploratory 30° navigation threshold.
7. Compare known experimental contrasts to explicitly labelled controls,
   keeping missing or inapplicable controls visible. No absent control is
   silently replaced with an invented same-state replicate.

## How to run

From the repository root, with Node.js, R and `jsonlite`:

```bash
node scripts/atlas-multiprotein-fetch.cjs benchmarks/output/atlas-multiprotein
Rscript scripts/atlas-multiprotein-benchmark.R benchmarks/output/atlas-multiprotein
```

The GitHub Actions workflow `Multi-protein experimental Atlas benchmark`
runs these commands on a fresh environment. The artifact includes source
manifest, coordinate records, per-residue data and a comparison report.

## First live-data results (2026-10-08)

The [live workflow](https://github.com/BiKC/RamplotR/actions/runs/37775734519)
successfully analyzed all seven updated-mmCIF files. The
[complete archived dataset and per-position CSVs](https://github.com/BiKC/RamplotR/actions/runs/37775734519/artifacts/11549363281)
are retained as an Actions artifact.

| Protein and compared chains | Comparison | Shared UniProt positions | Cα dRMSD (Å) | Complete φ/ψ pairs | ≥30° local changes |
| --- | --- | ---: | ---: | ---: | ---: |
| ADK, 4AKE A / 1AKE A | Open vs closed | 214 | 6.5077 | 212 | 37 |
| ADK, 4AKE A / B | Same crystal | 214 | 0.4022 | 212 | 16 |
| ADK, 1AKE A / B | Same crystal | 214 | 0.2513 | 212 | 10 |
| MBP, 1OMP A / 1ANF A | Open vs closed | 370 | 2.8060 | 368 | 23 |
| MBP, 1OMP A / 1JW4 A | Separate apo crystal | 370 | 0.4119 | 368 | 13 |
| RBP, 1URP A / 2DRI A | Open vs closed | 271 | 2.7443 | 269 | 10 |
| RBP, 1URP A / B | Same crystal | 271 | 0.0773 | 269 | 0 |

Three additional **self-comparisons**, one per protein, returned exactly zero
global dRMSD and zero local candidates. Each experimental comparison used the
same canonical residue positions for both structures. Every comparison had
zero mismatched deposited monomer identities among compared positions.
The synthetic cross-UniProt and substantial construct-mismatch cases were
rejected as intended.

**Method consistency:** Cα dRMSD is averaged over distinct residue pairs,
excluding the zero self-distance diagonal, matching the original ADK case
study. The observed-position overlap percentages are relative to the
two **observed, successfully mapped chains**, not the full UniProt sequence
or all crystallographically unresolved termini.

Three known contrasts are larger than the chosen matched controls by global
dRMSD, but the local angle-change fraction depends strongly on the protein:
ADK 37/212, MBP 23/368 and RBP 10/269. Same-state controls have 16/212
and 10/212 (ADK), 13/368 (MBP), and 0/269 (RBP) local changes. Consequently,
30° is useful as a **visual exploration threshold**, not a validated
biological-state boundary.

### Source provenance

PDBe updated-mmCIF source SHA256 values from the benchmark run:

| PDB | SHA256 |
| --- | --- |
| 4AKE | `3b1cad099974d4d78d218b5699480a453701533612e696d92d51b7621efdcebd` |
| 1AKE | `7ba43f1063e8f43efd138d79cff888ee427b5dc42b17ca5a35a46496c8d120e0` |
| 1OMP | `1f6631b739a038cbf3a8b0dbc9e65c9b66bed1396cafbe6e4fbc30bf4111c16e` |
| 1ANF | `97782731afcee202e64d0b6846a2344bc02eedb3cfbf3fa06a7ff6b7f86c4c7b` |
| 1JW4 | `5413313d3f3179614a8cb598e12fd57cb8e7c25a0e6912cc0f40da75eef55512` |
| 1URP | `0302dc530a4ad91f35ce3ffeb24cc117339fa67d1699ed74d5d1a4d1e593d8b4` |
| 2DRI | `cf3132729848d8eee53465f7334b1ae523f88c6c40cd833dd5223c7dc42cb6c1` |

## Limitations

Experimental conditions, ligands and crystal contacts are different. Even
a clear positive/negative difference in these curated examples does not imply
causation or general sensitivity and specificity. Same-crystal chain copies
are not independent replicates. Small local backbone-angle differences do
not imply no domain motion; global domain movement and backbone torsion are
separate observations.

The 1.5 Å exploratory Atlas geometry-group cutoff and 30° local angular
navigation threshold are **not tuned or evaluated as classifiers** here.
This benchmark will be expanded with independent folds, construct controls
and blind structure pairs before any calibration is justified.

## Research sources

- Adenylate kinase: <https://www.rcsb.org/structure/4AKE> and
  <https://www.rcsb.org/structure/1AKE>
- Maltose-binding protein: <https://pmc.ncbi.nlm.nih.gov/articles/PMC240646/>
  and <https://www.rcsb.org/structure/1JW4>
- Ribose-binding protein: <https://pmc.ncbi.nlm.nih.gov/articles/PMC8150535/>
  and <https://www.rcsb.org/structure/1URP>
