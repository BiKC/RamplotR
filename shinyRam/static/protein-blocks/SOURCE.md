# Protein Blocks reference provenance

RamplotR uses the 16-state Protein Blocks (PB) structural alphabet as a
fragment-level backbone representation for conformational comparison.

## Scientific definition

Protein Blocks were introduced by:

A. G. de Brevern, C. Etchebest and S. Hazout.
"Bayesian probabilistic approach for predicting backbone structures in terms
of protein blocks." Proteins 41:271-288 (2000).

Each PB represents a five-residue local backbone fragment using eight
dihedral angles in this order:

1. psi(n-2)
2. phi(n-1)
3. psi(n-1)
4. phi(n)
5. psi(n)
6. phi(n+1)
7. psi(n+1)
8. phi(n+2)

RamplotR assigns the PB with the smallest periodic angular RMSD. Missing
required angles or a peptide-chain break leave the PB unassigned.

## Reference implementation

The reference angle table and assignment convention were checked against the
public PBxplore implementation:

- project: https://github.com/pierrepo/PBxplore
- reference definitions: `pbxplore/PB.py`
- assignment: `pbxplore/assignment.py`
- documentation: https://pbxplore.readthedocs.io/

PBxplore wraps angular differences through +/-180 degrees before summing
squared deviations. RamplotR uses the same nearest-prototype criterion and
reports the square root of the mean squared wrapped angular difference as
`protein_block_rmsda`; taking the square root/mean does not change which
prototype is nearest.

The RamplotR implementation is written in R and uses the already calculated
RamplotR backbone torsions. It additionally requires peptide continuity across
all five residues when `bonded_to_next` is available.

## License notice for upstream reference material

PBxplore is distributed under the MIT License.

Copyright (c) 2013 Poulain, A. G. de Brevern

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the condition that the copyright and permission
notice are included in copies or substantial portions of the Software.

The upstream software is provided without warranty; consult the PBxplore
LICENSE file for the complete license text.

## Interpretation in RamplotR

Protein Blocks are used as a local conformational fingerprint. They are not:

- a stereochemical validation category;
- a DSSP/secondary-structure assignment;
- a statement about biological function;
- evidence of molecular dynamics or thermodynamic populations.

The simpler Alpha-R/Beta/PPII/Alpha-L/Other labels remain available as a
human-readable one-residue overview. Protein Blocks provide the higher
resolution fragment representation intended for Atlas switch-region analysis.
