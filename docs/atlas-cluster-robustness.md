# Experimental Atlas cluster sensitivity

Atlas groups experimental structures using exact PDBe SIFTS/UniProt-mapped
first-model C-alpha coordinates. A suggested geometric grouping can change if
part of the shared sequence is excluded. This check reports how sensitive an
existing grouping is to **which backbone positions are used**.

## How it works

For a selected cohort of at least three structures, Atlas:

1. Starts with the same common UniProt core used in the main distance-map
   RMSD and average-linkage dendrogram. Large cores are already subsampled
   deterministically to at most 300 ordered positions.
2. Divides the ordered sampled positions into ten consecutive blocks. When
   there are gaps in observed residues, a block follows the sampled sequence
   order, not an invented continuous sequence interval.
3. Omits each block once. For each omission, recalculates the entire cohort
   distance matrix from internal C-alpha pair distances on the remaining
   common positions.
4. Reruns the selected **automatic cluster suggestion** or the exact
   user-selected **manual distance cutoff**. In automatic mode, the candidate
   group count can change across omissions.
5. Compares partitions by their member pairs. Arbitrary numeric group names
   are ignored.

There is no random seed and no new sequence alignment or network download.
Repeating the same cohort with identical parameters reproduces the same
sensitivity table.

## Reading the results

- **Whole-partition recovery**: number of ten omissions that reproduce every
  original same/different-group relationship.
- **Exact group recovery**: number of omissions where a particular original
  group reappears with exactly the same members, independent of its label.
- **Pair co-assignment**: number/fraction of omissions where a specific pair
  of structures is grouped together. It can be high for pairs originally
  separated or low for originally co-grouped pairs.
- **Same PDB entry**: flags structures whose entity identifiers share the
  same four-character PDB deposition ID. These are not independent
  experiments just because they are separate chains or polymer entities.
- **Omitted-block runs**: start/end UniProt positions, omitted/retained count,
  resulting group count and agreement with the original partition.

All three tables are downloadable as CSV. The grouping itself and the cluster
representatives are unchanged by this diagnostic.

## Limitations

This is **not** a bootstrap confidence interval, sampling probability,
likelihood of an experimental state, or evidence of functional-state identity.
Adjacent canonical positions and multiple polymer entities from the same
deposition are not independent biological replicates. The check cannot test
unobserved regions, missing conformational states, ligand annotations,
experimental method biases or construct equivalence. A stable single-cluster
result may simply reflect conservative automatic cut safeguards. For a cohort
of two structures, the sensitivity score is deliberately unavailable.

Construct and monomer-chemistry compatibility is separately assessed before
clustering and still requires explicit acknowledgement when chemistry differs.
Do not apply arbitrary robustness thresholds to declare functional states.

## Implementation and testing

- Algorithm: `shinyRam/R/atlas-robustness.R`
- Integration: `ram_atlas_geometry_groups()` with its current verified
  canonical-core matrices; both manual and automatic modes are supported.
- Deterministic, shared-deposition, no-split and invalid-input controls:
  `tests/atlas-geometry.R`
- Browser notice for a two-entry cohort: `tests/atlas-browser.cjs`

The current algorithm is intended for interactive cohorts of 2–12 verified
polymer entities, not a large archive-wide bootstrap study.
