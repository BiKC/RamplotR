# Experimental Atlas benchmark: adenylate kinase

This case study compares experimental structures of *Escherichia coli*
adenylate kinase, not predicted conformational populations.

## Inputs and provenance

- **4AKE**: open, apo, crystallographic adenylate kinase.
- **1AKE**: closed, Ap5A inhibitor-bound, crystallographic adenylate kinase.
- Both refer to UniProt **P69441**. Exact residue mapping is extracted directly
  from PDBe's updated mmCIF (`_pdbx_sifts_xref_db` and
  `_pdbx_poly_seq_scheme`).
- The same `atlas-sifts-exact.js` parser used by the browser extracts SIFTS
  correspondences and first-model backbone coordinates.
- Do not assume an author residue number is the UniProt position.
- Every live run writes the downloaded mmCIF SHA256 checksum and source URL
  into a machine-readable `manifest.json`.

## Comparisons

1. 4AKE's first covered chain against 1AKE's first covered chain (open/closed).
2. Two chain copies within 4AKE (open-state crystal packing/chain control).
3. Two chain copies within 1AKE (closed-state crystal packing/chain control).
4. A chain compared to itself (numerical zero control).

These chain comparisons are **not independent experiments or biological
replicates**. A single known open/closed pair cannot estimate classification
sensitivity, specificity or generalizable cutoff values.

## Readouts

- Strict common canonical UniProt Cα coverage: at least 100 positions and
  at least 60% of both chains' mapped positions.
- Intrachain Cα distance-map RMSD, invariant to rigid orientation.
- Per-UniProt-position differences in all-to-all Cα distances, an indicator
  of global/domain rearrangement.
- Circular φ/ψ differences only where both chains have continuous peptide
  geometry and exact canonical correspondence.
- Number of comparable φ/ψ positions, 30° navigation-threshold candidates,
  and contiguous candidate switch regions.
- Descriptive CORE, LID and NMP means for *distance-map contributions*.
  These regions are approximate literature-based labels for this example,
  not a structural-state ground truth for threshold training.

## Run

From the repository root with Node.js, R and `jsonlite` installed:

```bash
node scripts/atlas-adk-fetch.cjs benchmarks/output/atlas-adk
Rscript scripts/atlas-adk-benchmark.R benchmarks/output/atlas-adk
```

Alternatively, run the **Experimental Atlas biological benchmark** workflow
in GitHub Actions. It uploads the exact residue tracks, comparison summaries,
source manifest and Markdown report as an artifact.

The live benchmark is intentionally separate from the regular offline unit
tests. PDBe outages or schema changes should fail the live benchmark clearly
without blocking unrelated offline scientific regression tests.

## First live-data results (2026-10-08)

GitHub Actions run [37765740664](https://github.com/BiKC/RamplotR/actions/runs/37765740664)
retrieved the PDBe **updated mmCIF** files and completed both structure
comparisons. The [benchmark artifact](https://github.com/BiKC/RamplotR/actions/runs/37765740664/artifacts/11544028276)
contains the exact source manifest, atom/mapping records, tracks and reports.

| Compared model-1 chains | Common UniProt positions | Cα distance-map RMSD (Å) | Complete φ/ψ pairs | Above 30° | Candidate segments |
| --- | ---: | ---: | ---: | ---: | ---: |
| 4AKE A vs 1AKE A (open/closed) | 214 | 6.5077 | 212 | 37 | 22 |
| 4AKE A vs 4AKE B (same crystal) | 214 | 0.4022 | 212 | 16 | 9 |
| 1AKE A vs 1AKE B (same crystal) | 214 | 0.2513 | 212 | 10 | 8 |
| 4AKE A vs itself | 214 | 0.0000 | 212 | 0 | 0 |

Mean per-residue contributions to the *global Cα distance-map change*
between 4AKE A and 1AKE A (Å): CORE 4.540, LID 8.600, NMP 7.591.
This is distinct from the local circular φ/ψ metric.

Both structures mapped to 214 observed UniProt positions across two chains,
with 428 Cα and 1,284 complete N/Cα/C atom records per entry.

Input source SHA256 digests from that run:

- `4ake_updated.cif`:
  `3b1cad099974d4d78d218b5699480a453701533612e696d92d51b7621efdcebd`
- `1ake_updated.cif`:
  `7ba43f1063e8f43efd138d79cff888ee427b5dc42b17ca5a35a46496c8d120e0`

The global open/closed difference clearly exceeds the two within-crystal
controls, as expected for this documented domain-motion pair. However,
the chain-copy controls still have 16 and 10 residues above the 30°
navigation threshold. A fixed angular cutoff therefore cannot be
interpreted as a validated positive/negative functional-state classifier.

The current benchmark check protects against a **loss of contrast on this
specific pair**; it is not a trained, transferable classification threshold.

## Interpretation limits

Open/closed adenylate kinase is a documented large-domain-motion example.
A large global geometry difference can coexist with modest local φ/ψ
changes over most of the sequence. A locally changed backbone does not by
itself establish catalytic mechanism, a ligand-caused change, or an error in
either experimental structure.

This case study checks the mechanics of RamplotR on biological data. It is
**not yet** a validated protein conformational-state classifier. Future work
requires multiple proteins, independent experiments, construct compatibility
reviews, same-state positive controls and error analyses.

## Literature

- Müller et al., *Structure*, 1996. Experimental comparison of unligated
  *E. coli* adenylate kinase with an inhibitor-bound structure.
- Beckstein et al., *J Mol Biol*, 2009. Adenylate kinase open/closed
  conformational transitions (4AKE and 1AKE).
- Current experimental records: https://www.rcsb.org/structure/4AKE and
  https://www.rcsb.org/structure/1AKE.
