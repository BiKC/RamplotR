# Canonical UniProt residue coordinates

RamplotR keeps the original PDB author chain, residue number and insertion code
as the local coordinate system used for structure inspection. It can also add a
canonical UniProt coordinate when that mapping is available.

This second coordinate system is the foundation for the Conformational Atlas:
structures with different constructs or residue numbering need a shared protein
coordinate system before they can be aggregated safely.

## Public PDB accessions

For structures loaded by PDB accession, RamplotR requests the PDBe SIFTS
PDB-to-UniProt mapping in the browser. PDBe supports browser/AJAX access, so the
same request path works in the normal Shiny application and in Shinylive without
adding a server-side HTTP dependency.

The mapping response contains UniProt accessions and chain mapping segments with
PDB/author and UniProt range endpoints, sequence identity and coverage.

RamplotR stores those segment records separately from the residue table.

### Conservative range expansion

A SIFTS segment is expanded into residue-level canonical coordinates only when
all of the following are true:

- author start/end residue numbers are available;
- UniProt start/end positions are available;
- neither range endpoint has an insertion code;
- both ranges proceed forward;
- the author-number span and UniProt span are exactly equal.

For example:

    PDB chain A author residues 101-103
    UniProt P12345 residues 10-12

can be expanded safely as 101→10, 102→11, 103→12.

A range such as 200-201A, or ranges with unequal spans, is retained as a SIFTS
segment but is **not** converted into invented residue positions.

This is intentionally conservative. A later Atlas phase will add exact
residue-level SIFTS data from PDBe-enriched mmCIF/SIFTS sources for nonlinear
segments.

## AlphaFold DB

AlphaFold DB entries are requested by UniProt accession and use the model's
UniProt sequence numbering. RamplotR therefore adds that accession and residue
position directly to the single protein chain.

This direct mapping applies specifically to AlphaFold DB input. RamplotR does
not assume that an arbitrary uploaded predicted structure uses canonical
UniProt numbering.

## User interface

When a selected residue has an unambiguous canonical mapping, the shared
inspector shows both coordinates, for example:

    Chain A 101 · SER
    UniProt P12345:10

The residue table also adds UniProt accession/position columns when mappings are
available. The Summary view reports how many currently selected residues are
mapped.

Ambiguous mappings are marked as ambiguous rather than assigned a guessed
canonical position.

## Provenance and limitations

Canonical coordinates are annotations. They do not replace:

- PDB author residue numbers;
- insertion codes;
- chain identifiers;
- sequence-alignment evidence used by pairwise comparison.

SIFTS is maintained by PDBe/UniProt and is the intended common residue mapping
for cross-structure aggregation. Network retrieval can fail or an entry may
contain a mapping pattern that cannot yet be expanded safely. In those cases
RamplotR continues to work with local PDB coordinates and reports the canonical
mapping as unavailable or partial.

The Conformational Atlas should use canonical positions only where mapping
provenance and coverage are explicit.
