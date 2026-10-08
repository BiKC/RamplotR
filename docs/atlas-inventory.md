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
cross-reference rows. Ambiguous and unresolved positions are not silently
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
- State clustering and representative-state selection are not yet computed.
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
