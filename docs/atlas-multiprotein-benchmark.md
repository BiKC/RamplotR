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
