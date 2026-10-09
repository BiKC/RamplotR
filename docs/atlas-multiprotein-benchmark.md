# Atlas experimental benchmark across four proteins and three fold families

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
| Pig-heart citrate synthase (P00889) | 1CTS, open | 2CTS, product-bound closed | 1CTS vs 3ENJ, separate crystal open-like control; 3ENJ contains cystamine-related covalent chemistry |

One exact self-comparison is added for every protein. It must have zero
geometry difference and zero local-angle changes.

Maltose-binding and ribose-binding protein are distinct proteins but belong
to related periplasmic binding-protein structural families. The citrate
synthase case adds a distinct all-alpha enzyme fold in addition to the
kinase and periplasmic-binding families. Three fold-family groupings are
still not a representative sample of structural diversity.

The pig-heart citrate synthase case is supported by the
[original open/closed structures](https://pdb101.rcsb.org/motm/93) and
[3ENJ experimental paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC2675578/).
This paper specifically identifies 1CTS and 3ENJ as open-like, and 2CTS
as closed, while documenting the covalently modified Cys184 in 3ENJ.
That modification and differences in experimental conditions are confounders,
not a claim of identical constructs or ligand states.

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

## Expanded fourth-fold benchmark (2026-10-09)

The curated manifest now includes 1CTS, 2CTS and 3ENJ and labels the fold
family for all four proteins. This raises the benchmark to **ten experimental
PDB entities, four proteins and three structural fold-family groupings**.

The live workflow reports an additional `per-protein-control-contrasts.csv`
with the documented contrast, strongest available measured control,
separate-crystal versus same-crystal control counts, and both global
distance-map and local φ/ψ difference margins. These are **descriptive
measurements, not threshold-based accuracy or statistical validation**.

The first three-protein run remains an immutable historical snapshot below.
The fourth-family extension completed its
[live PDBe CI run](https://github.com/BiKC/RamplotR/actions/runs/37909618876)
on October 9, 2026. The immutable run artifact contains exact SIFTS
records, model-1 coordinates, per-residue torsion/geometry readouts and
full PDB source-file SHA256 values.

| Citrate synthase comparison | Shared exact UniProt positions | Cα dRMSD (Å) | Complete φ/ψ pairs | ≥30° candidates | Deposited monomer-identity differences |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1CTS A (open) / 2CTS A (closed) | 437 | 1.6241 | 435 | 91 | 1/437 |
| 1CTS A (open) / 3ENJ A (independent open-like) | 437 | 0.6962 | 435 | 103 | 0/437 |

**Important counterexample:** despite the greater global dRMSD in the
documented open/closed comparison, the independently crystallized open-like
control shows **more** ≥30° local angular changes (103 versus 91).
This directly contradicts treating the local angle threshold as a
general functional-state discriminator. The one differing experimental
monomer in 1CTS/2CTS and the covalent Cys184 chemistry of 3ENJ must not
be treated as proven irrelevant to the conformational observations.

The report intentionally records a *descriptive* global distance-map
contrast margin (1.6241 − 0.6962 = 0.9279 Å) and a negative local
candidate-fraction contrast (91/435 versus 103/435). These are not effect
sizes with inferential uncertainty or generalizable classifier outcomes.

**New updated-mmCIF SHA256 provenance:**

| PDB | Source SHA256 |
| --- | --- |
| 1CTS | `10886ce9eb0597314aa9211c0771dd0a9a20eb9c848cbf5a28bfd42bb0160ab3` |
| 2CTS | `a4620d15d23b22554a6fd42a9ff395c8df0e8de54f00555feb37f08a843e9120` |
| 3ENJ | `5265b3083f1c765f74fab80531d6687a8fcafd8408fc7660ab88795edc985a02` |

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

In the original three-protein set, three known contrasts were larger than the chosen matched controls by global
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
Additional independent-fold cases, construct/ligand controls and blinded
structure pairs are still required before any calibration is justified.
One new fold and one independent open-like control do not establish
population-level sensitivity, specificity or a universal state threshold.

## Research sources

- Adenylate kinase: <https://www.rcsb.org/structure/4AKE> and
  <https://www.rcsb.org/structure/1AKE>
- Maltose-binding protein: <https://pmc.ncbi.nlm.nih.gov/articles/PMC240646/>
  and <https://www.rcsb.org/structure/1JW4>
- Ribose-binding protein: <https://pmc.ncbi.nlm.nih.gov/articles/PMC8150535/>
  and <https://www.rcsb.org/structure/1URP>
- Citrate synthase: <https://pmc.ncbi.nlm.nih.gov/articles/PMC2675578/>,
  <https://pdb101.rcsb.org/motm/93>, and
  <https://www.rcsb.org/structure/2CTS>
