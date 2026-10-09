# Experimental Conformational Atlas inventory

The **Atlas** tab currently performs experimental structure discovery for a
UniProt accession. It is the first retrieval step towards a full conformational
atlas, not yet a classifier of structural states.

## Use

1. Open the Atlas tab.
2. Enter a UniProt accession, for example `P00533`.
3. Choose **Find experimental structures**.
4. Inspect the returned experimental PDB polymer entities, chains, method,
   release date, resolution and reported UniProt sequence coverage.
5. With a primary structure loaded, select **Compare with loaded structure**
   to inspect a candidate in the existing Conformational Change Explorer.

The search uses an exact UniProt cross-reference attribute in the RCSB PDB
Search API and requests experimental results only. For each returned entity,
the browser reads public entity and entry metadata from the RCSB Data API.
The first response is capped at 50 polymer entities; the displayed total hit
count distinguishes a partial inventory from a complete result. At most six entity-enrichment tasks run concurrently (each task may fetch
both entity and entry metadata). Failed metadata fetches are counted,
rather than silently presented as complete records.

Click **Verify SIFTS mapping** on one experimental entity to download the
PDBe **updated mmCIF** and validate its residue-level UniProt correspondence.
This is an explicit action because mmCIF files can be large. RamplotR joins
the SIFTS `asym_id/seq_id` pair to the polymer sequence scheme, then preserves
PDB chain ID, residue number, insertion code and UniProt accession/position.
Results report exact mapped positions, observed residues, conflicts and unmatched
cross-reference rows. Once at least one entity has been verified, the Atlas
shows an experimental cohort summary and two CSV exports:

- exact residue rows with PDB entry, entity, author/PDB chain/number/insertion,
  UniProt accession/number, observed state, provenance and ambiguity status;
- per-UniProt-position support, counting distinct **verified PDB entities**
  with an observed, unambiguous mapping.

These counts are **not** the number of independent experimental measurements,
protein conformations, or sequence coverage.

## Experimental geometry groups (exploratory)

When two or more verified experimental entities each have at least 30
unambiguous, observed C-alpha positions, the Atlas offers
**Compare structural geometries**. It uses the first coordinate model in
the downloaded PDBe updated mmCIF. Each coordinate must match an exact
SIFTS polymer `label_asym_id/label_seq_id` reference. Alternate atom
conformers, missing atoms, ambiguous mappings and secondary protein copies
are not counted as extra observations.

The comparison requires at least 30 shared observed UniProt residues and
at least 60% common-core coverage of each selected chain. The common
position set is identical for all compared entities. For performance,
large common cores are sampled evenly to at most 300 positions; the number
of positions used is always displayed.

RamplotR calculates the root mean squared difference between intrachain
C-alpha distance maps (distance-map RMSD, Å). This captures internal
geometry changes without needing a rigid-body superposition. Average-linkage hierarchical clustering provides **exploratory geometric groups**,
with one medoid per suggested group. Automatic clustering is the default and
may suggest a single group when the available structures do not support a
clear split. The manual distance cutoff remains available (initially 1.5 Å).
Neither mode assigns functional-state labels.

These groups must not be presented as functional states, experimental
populations, independent structural observations, ligand-caused transitions,
or calibrated conformational classes. Different constructs and conformations
can affect distances, and a global distance map does not identify which
local backbone residues change. Those checks are separate roadmap steps.

The geometry analysis uses data from explicitly verified structures only.
It does not download every experimental entry from the RCSB inventory.

## Local backbone-change candidates

After comparing experimental geometry groups, RamplotR can compare the first
two group representatives at their exact UniProt-mapped positions. It extracts
model-1 **N, Cα and C** backbone atoms from the same PDBe updated mmCIF file.
It calculates φ/ψ only when consecutive native polymer sequence positions
have complete backbone atoms and a 1.0–1.9 Å C–N peptide connection. Gaps,
ambiguous atoms and mapping conflicts never generate extrapolated torsions.

Circular Δφ/Δψ values are wrapped independently across ±180°. A combined
angular displacement is used solely to find noteworthy local differences.
The default 30° navigation threshold can be changed in the Atlas. Consecutive
canonical positions above threshold are grouped into *candidate change
segments*, with isolated residues reported separately as one-position
segments. The plot and residue CSV retain missing-comparison positions.

The local inspection can compare any two explicitly selected verified experimental structures, including when Atlas suggests only one geometric cluster. The group representatives are used as defaults when more than one cluster exists. This local inspection still compares only two structures,
**not** all structures within each group, and does not perform a statistical
group test. A difference in φ/ψ does not prove a biological functional
transition, ligand effect, protein-block assignment, or error in either
experimental model. A rigid domain shift can also occur with little or no
local torsion change. Neither the geometric distance cutoff nor the local
30° threshold has yet been calibrated on curated biological case studies. Ambiguous local residue mappings
remain in the residue export but never count as canonical position support. Ambiguous and unresolved positions are not silently
assigned a canonical coordinate.

The RCSB UniProt-reference-coverage field is entity-level metadata, **not**
the verified per-residue canonical mapping. It does not prove two structures
represent the same protein state. Coverage unavailable from the API is shown
as unavailable.

## Limitations

- One polymer entity may occur in multiple chains within an entry.
- Separate experimental entries may reproduce essentially the same structural
  state.
- Differences in constructs, sequence variants and bound partners require
  careful interpretation.
- Results can be paginated through the currently reported archive hits, but
  do not represent a frozen archive snapshot. If the reported hit count
  changes mid-search, restart rather than merging potentially incompatible
  pages.
- Page offsets count raw RCSB search hits. Metadata failures and duplicate
  entity IDs are reported separately; neither silently shifts pagination.
- When entity metadata is unavailable, a disclosure lists the affected
  PDB/entity identifiers so researchers can distinguish a failed enrichment
  from an absent experimental structure.
- Geometric clustering and medoid selection are implemented, but these are not validated biological-state assignments.
- Exact residue mapping is now available on demand for selected entities.
  Full cohort-scale mapping, isoform equivalence and consistent coverage
  filtering before state clustering remain future work.
- PDBe updated mmCIF uses `_pdbx_sifts_xref_db` and
  `_pdbx_poly_seq_scheme`; the browser must be able to reach PDBe.
  Retrieval failures leave the inventory intact.
- Selecting a structure for comparison is a user-controlled network action.
  Arbitrary uploaded structures are not transmitted to the search service.

The eventual Conformational Atlas will add canonical coverage checks,
structural-state clustering, representative structures and local conformational
switch regions. See [the roadmap](conformational-atlas-roadmap.md).

## Regression tests

- `node tests/atlas-discovery.test.cjs` tests normalized metadata and bounded
  enrichment concurrency without network access.
- `node tests/atlas-browser.cjs` exercises a live local Shiny application
  against mocked, multi-page RCSB search/metadata responses and captures
  desktop/mobile screenshots. It runs in the browser-preview GitHub Actions
  workflow.
- `Rscript tests/atlas.R` checks page order, archive-total changes,
  deduplication and failure recovery in the pure-R cohort merge helper.
