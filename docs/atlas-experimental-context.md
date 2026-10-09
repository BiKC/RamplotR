# Experimental conditions and confounders in Atlas

Experimental structures can be geometrically similar but differ in sample
preparation, experimental method, deposition provenance or amino-acid chemistry.
The current Atlas cannot infer a biological state or a ligand-induced mechanism
from geometry.

## Available provenance

The same RCSB Polymer Entity and Entry responses already retrieved by the
experimental search supply:

- experimental method(s);
- reported structure resolution, when available;
- initial PDB release date;
- polymer-entity description.

The Atlas geometry step already verifies observed UniProt position identity and
audits deposited polymer monomer chemistry through exact PDBe updated-mmCIF
SIFTS mappings. The experimental-context panel joins these two sources **by
PDB code plus polymer entity number** and associates each selected entity with
its exploratory geometric group.

## Reading the context audit

The per-entity table shows method, resolution, release date, observed canonical
coverage and fraction of observed positions with known deposited monomer IDs.
Unavailable metadata stays missing.

The pairwise table reports:
- whether two structures appear in the same Atlas geometric group;
- whether both polymer entities come from one PDB deposition (not independent
  experimental replicates);
- differences in listed experimental methods;
- differences in resolution when both values exist;
- documented discrepancies and unknowns in polymer monomer chemistry;
- differences in observed sequence spans.

A PDB release date is **not** a measurement of when a protein conformation
occurred. Resolution is reported, not used to rank biological states, and it
is not directly comparable across experimental methods without context.

Importantly, **ligand/cofactor occupancy is not verified by this workflow**.
The panel does not call a group apo, holo, WT, mutant, active, inactive or
experimentally independent. Polymer-chemistry differences can also reflect
modified residues rather than sequence mutations.

## Exports and validation

Both selected-entity and pairwise metadata can be exported as CSV. Joining
requires the same UniProt cohort and exact selected PDB/entity identifiers.
Missing, duplicate and mismatched archive records are not guessed or borrowed
from other structures.

Source code: `shinyRam/R/atlas-experimental-context.R`.
Standalone tests: `tests/atlas-experimental-context.R`.
