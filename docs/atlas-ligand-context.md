# Observed non-water components near experimental structures

The Atlas can report components recorded as non-water HETATM residues in
the same PDBe updated mmCIF file used for exact SIFTS/UniProt verification.
The aim is to help investigators check ligand/cofactor *context* without
automatically labelling a PDB entry apo, holo or functionally active.

## Measurement

For every verified polymer entity in a clustered Atlas cohort:

1. Read model 1's atom_site heavy-atom coordinates from the already
   downloaded updated mmCIF file, using exact PDB label asym and sequence
   identifiers to select the observed UniProt-mapped target chain.
2. Exclude water, hydrogen and deuterium, non-first models, and alternate
   conformers other than blank or A. Modified polymer residues with defined
   label sequence positions are **not** treated as free ligands.
3. Identify non-water, nonpolymer HETATM residues and find those with a heavy
   atom within 4.5 Å of any mapped heavy protein atom. Only positions with
   unambiguous exact mapping are used.
4. Report deposited three-letter component code, ligand asym/author residue
   identifiers, nearest exact UniProt residue and minimum heavy-atom distance.

The spatial search uses 5 Å bins and inspects neighboring bins. It never
approximates contacts from residue numbering or uses AlphaFold pLDDT as a
substitute. The method is a proximity search, not a chemical binding model.

## Interpretation

- **Measured and nearby** means a deposited non-water HETATM component has
  at least one heavy atom within the cutoff of the selected target chain.
  This can include ions, crystallization additives and nonfunctional contacts.
- **Measured, no nearby site** means the extraction completed but no such
  deposit record was within the cutoff. It does not prove an apo structure.
- **Unavailable** means extraction failed or the file lacks appropriate
  coordinates; the structure may still have valid SIFTS verification and
  contribute to geometric clustering.
- **No recorded non-water HETATM site** is a property of the examined
  model's deposited coordinates, not experimental confirmation of ligand
  absence.

Occupancy, biochemical identity, covalent chemistry, protonation, assembly,
protein-binding partners and genuine ligand-dependent functional changes
are not established. Ligands represented as polymer components, unmodelled
ligands, excluded alternates and ligands near unmapped residues can be missed.
The result is specific to model 1 and the selected exact-mapped chain.

Researchers should inspect the structure paper, experimental conditions and
wwPDB chemical component description before calling an observed compound
an activating, inhibitory or physiologically bound ligand.

## Files and reproducibility

- `shinyRam/www/atlas-ligand-contacts.js`: targeted atom_site parser and
  nearest mapped heavy-atom calculation.
- `shinyRam/R/atlas-ligand-evidence.R`: provenance-labelled cohort table.
- `tests/atlas-ligand-contacts.test.cjs`: contact, solvent, modification,
  model, alternate location and missing-coordinate controls.
- `tests/atlas-ligand-evidence.R`: strict separation of measured,
  no-local-contact and unavailable.
- `tests/atlas-browser.cjs`: end-to-end synthetic experimental PDB
  upload/download verification and Atlas panel.

Two CSV exports record contact availability for every selected structure
and individual deposited components meeting the fixed 4.5 Å criterion.
No new third-party network endpoint is necessary.
