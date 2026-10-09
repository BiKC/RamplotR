# Standalone RCSB metadata + verified-geometry confounder regression.
source(file.path("shinyRam","R","canonical.R"))
source(file.path("shinyRam","R","atlas.R"))
source(file.path("shinyRam","R","atlas-experimental-context.R"))
assert <- function(x,message) if(!isTRUE(x))stop(message,call.=FALSE)
ids <- c("1ABC_1","1ABC_2","2XYZ_1")
cohort <- list(accession="P12345",results=list(
  list(pdb_id="1ABC",entity_id="1",method="X-RAY DIFFRACTION",
    resolution=1.9,release_date="2021-01-12",description="Test protein"),
  list(pdb_id="1ABC",entity_id="2",method="X-RAY DIFFRACTION",
    resolution=1.9,release_date="2021-01-12",description="Test protein"),
  list(pdb_id="2XYZ",entity_id="1",method="ELECTRON MICROSCOPY",
    resolution=3.2,release_date="2023-05-17",description="Test protein")
))
geometry <- list(accession="P12345",selected=ids,assignment=data.frame(
  entity=ids,geometric_group=c(1L,1L,2L),
  common_coverage=c(.95,.95,.9)))
pairs <- data.frame(entity_a=c(ids[1],ids[1],ids[2]),
  entity_b=c(ids[2],ids[3],ids[3]),common_observed=c(40,40,40),
  chemistry_differences=c(0,1,0),chemistry_unknown=c(0,0,5),
  unmatched_a=c(0,0,0),unmatched_b=c(0,4,4))
profiles <- stats::setNames(replicate(3L,
  list(known=40L,observed=40L),simplify=FALSE),ids)
audit <- list(selected=ids,pairs=pairs,profiles=profiles)
out <- ram_atlas_experimental_context(cohort,geometry,audit)
assert(nrow(out$entries)==3L && nrow(out$pairs)==3L &&
       out$shared_pdb_pairs==1L &&
       out$differing_method_pairs==2L &&
       out$differing_chemistry_pairs==1L,
  "Experimental context should report method/chemistry differences and shared PDB.")
assert(out$pairs$same_geometric_group[[1L]] &&
       !out$pairs$same_geometric_group[[2L]] &&
       out$pairs$same_pdb_entry[[1L]] &&
       grepl("Same PDB deposition",out$pairs$context_warning[[1L]]),
  "Same-deposition and cross-cluster comparison must remain distinct concepts.")
assert(isTRUE(all.equal(out$pairs$resolution_gap_A[[2L]],1.3)),
  "Known resolutions should be compared without treating differences as scores.")
assert(grepl("no apo/holo assignment",out$ligand_status,fixed=TRUE),
  "Absence of ligand verification must be explicit.")
# Partial/missing RCSB enrichment cannot become invented experimental metadata.
cohort_missing <- cohort
cohort_missing$results <- cohort_missing$results[1:2]
missing <- ram_atlas_experimental_context(cohort_missing,geometry,audit)
assert(missing$missing_metadata==1L &&
       is.na(missing$entries$method[[3L]]) &&
       is.na(missing$entries$resolution_A[[3L]]) &&
       missing$unresolved_method_pairs==2L,
  "Missing enrichment should propagate as unknown, not assume a method.")
# Name/order of metadata is irrelevant; join by validated entry IDs.
reordered <- cohort
reordered$results <- rev(reordered$results)
permuted <- ram_atlas_experimental_context(reordered,geometry,audit)
assert(identical(permuted$entries,out$entries) &&
       identical(permuted$pairs,out$pairs),
  "Metadata must join by PDB+entity ID, never page order.")
wrong_cohort <- cohort
wrong_cohort$accession <- "Q99999"
assert(inherits(try(ram_atlas_experimental_context(
  wrong_cohort,geometry,audit),silent=TRUE),"try-error"),
  "Mixed accessions must not produce a confounder report.")
duplicated <- cohort
duplicated$results[[4L]] <- duplicated$results[[1L]]
assert(inherits(try(ram_atlas_experimental_context(
  duplicated,geometry,audit),silent=TRUE),"try-error"),
  "Ambiguous repeated archive metadata must be rejected.")
message("Atlas experimental context and confounder tests passed.")
