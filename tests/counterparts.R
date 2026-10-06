source(file.path("shinyRam","R","counterparts.R"))
assert <- function(x,msg) if(!isTRUE(x)) stop(msg,call.=FALSE)

assert(identical(ram_uniprot_accession(" p69905 "),"P69905"),
       "UniProt accessions must be normalized.")
assert(grepl("/uniprot/summary/P69905.json$",
             ram_3dbeacons_summary_url("P69905")),
       "3D-Beacons summary endpoint changed unexpectedly.")
assert(inherits(try(ram_uniprot_accession("not an id"),silent=TRUE),"try-error"),
       "Invalid UniProt accessions must be rejected.")

mock <- list(
  uniprot_entry=list(ac="P69905",sequence_length=142),
  structures=list(
    list(summary=list(
      model_identifier="1abc.1.A",model_category="EXPERIMENTALLY DETERMINED",
      provider="PDBe",coverage=0.91,resolution=1.8,
      uniprot_start=2,uniprot_end=131,experimental_method="X-RAY DIFFRACTION",
      model_url="https://example.org/1abc.cif",
      entities=list(list(entity_type="POLYMER",chain_ids=list("A"))))),
    list(summary=list(
      model_identifier="2xyz.1.B",model_category="EXPERIMENTALLY-DETERMINED",
      provider="pdbe",coverage=0.72,resolution=NA,
      uniprot_start=10,uniprot_end=110,
      entities=list(list(entity_type="POLYMER",chain_ids=list("B","C"))))),
    # Duplicate PDB/chain with worse coverage should be removed.
    list(summary=list(
      model_identifier="1abc.2.A",model_category="EXPERIMENTALLY DETERMINED",
      provider="PDBe",coverage=0.50,resolution=2.6,
      entities=list(list(entity_type="POLYMER",chain_ids=list("A"))))),
    # Predicted/non-PDBe records must never be presented as experimental.
    list(summary=list(
      model_identifier="AF-P69905-F1",model_category="AB-INITIO",
      provider="AlphaFold DB",coverage=1)),
    list(summary=list(
      model_identifier="3def.1.A",model_category="EXPERIMENTALLY DETERMINED",
      provider="OtherProvider",coverage=1))
  )
)

parsed <- ram_parse_experimental_counterparts(mock,"p69905")
assert(nrow(parsed)==2L,
       "Only unique experimentally determined PDBe records should remain.")
assert(identical(parsed$pdb_id,c("1ABC","2XYZ")),
       "Counterparts should rank by coverage and parse PDB IDs.")
assert(parsed$coverage[[1L]]==0.91 && parsed$resolution[[1L]]==1.8,
       "Counterpart coverage/resolution parsing changed.")
assert(identical(parsed$chains,c("A","B, C")),
       "Polymer chain labels were not preserved.")
assert(parsed$uniprot_start[[1L]]==2L && parsed$uniprot_end[[1L]]==131L,
       "UniProt coverage range changed.")

seen <- NULL
result <- ram_lookup_experimental_counterparts(
  "P69905",max_results=1L,
  fetch=function(url) { seen <<- url; mock })
assert(nrow(result)==1L && result$pdb_id[[1L]]=="1ABC" &&
       grepl("P69905.json",seen,fixed=TRUE),
       "Lookup must use the normalized accession and result limit.")

empty <- ram_parse_experimental_counterparts(list(structures=list()),"P69905")
assert(nrow(empty)==0L && "pdb_id" %in% names(empty),
       "No-hit responses must return a typed empty table.")

message("Experimental counterpart lookup tests passed.")
