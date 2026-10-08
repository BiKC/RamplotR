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
