# Experimental counterparts for predicted structures

RamplotR can discover experimentally determined PDB structures for a predicted
protein and load a selected structure directly into the pairwise comparison
workflow.

## Data source

Counterpart discovery uses the 3D-Beacons Hub API with a UniProt accession.
RamplotR only presents records whose provider is **PDBe** and whose model
category is explicitly **EXPERIMENTALLY DETERMINED**.

This distinction is deliberate: template-based models, AlphaFold models,
conformational ensembles and records from other providers are not relabelled as
experimental evidence.

## Workflow

For an AlphaFold DB model, RamplotR already knows the UniProt accession and
performs the lookup automatically after the model is loaded.

For an uploaded AlphaFold 2 / ColabFold, ESMFold or other predicted structure,
open the prediction-confidence panel, enter the matching UniProt accession and
choose **Find experimental structures**.

Returned structures are ranked by:

1. mapped UniProt sequence coverage, highest first;
2. experimental resolution, lowest first when available.

Resolution can legitimately be unavailable for some experimental methods, such
as NMR; these structures remain in the result set.

Select one hit and choose **Compare selected structure**. RamplotR loads the PDB
entry into the existing Compare tab, where chain selection, sequence alignment,
wrapped delta-phi/delta-psi, the Conformational Change Explorer and linked
paired 3D inspection remain available.

## Interpretation

A structure mapped to the same UniProt accession is a useful experimental
counterpart, but it is not automatically an equivalent biological state.
Deposited structures may differ in:

- construct boundaries;
- mutations or engineered residues;
- ligand or cofactor state;
- oligomeric state;
- experimental conditions;
- unresolved segments;
- post-translational modifications;
- conformational state.

RamplotR therefore uses the lookup as a discovery/navigation step. Scientific
interpretation should still inspect sequence alignment, coverage, structure
metadata and the local residue context.

## Privacy and networking

Local uploaded coordinates are **not** sent to 3D-Beacons. The lookup sends only
the UniProt accession supplied by the user. Choosing a counterpart then loads
that public PDB entry through RamplotR's existing PDB input workflow.

## Reproducibility

The parsing and filtering logic is covered by `tests/counterparts.R` using a
mock 3D-Beacons response. Browser workflows do not depend on the live service,
so external downtime does not make RamplotR's core CI flaky.
