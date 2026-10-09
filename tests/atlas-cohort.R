source(file.path("shinyRam","R","canonical.R"))
source(file.path("shinyRam","R","atlas-sifts.R"))
source(file.path("shinyRam","R","atlas-cohort.R"))
assert <- function(x,msg) if(!isTRUE(x)) stop(msg,call.=FALSE)

make <- function(entity,positions) {
  map <- data.frame(chain="A",resi=seq_along(positions)+100L,
    insertion_code="",uniprot_accession="P12345",
    uniprot_resi=as.integer(positions),entity_id=entity,
    struct_asym_id="X",identity=NA_real_,coverage=NA_real_,
    observed=TRUE,canonical_source="PDBe updated mmCIF exact SIFTS")
  list(state="mapped",mapping=map)
}
verified <- list("1ABC_1"=make(1L,c(40L,41L,43L)),
                 "2XYZ_1"=make(1L,c(41L,43L,44L)),
                 "3ABC_1"=list(state="error",message="Unavailable"))
cohort <- ram_atlas_cohort_table(verified,"P12345")
assert("mon_id" %in% names(cohort) && all(is.na(cohort$mon_id)),
  "Legacy exact mapping without polymer chemistry must export unknown, not invent an identity.")
verified_with_chemistry <- verified
verified_with_chemistry[["1ABC_1"]]$mapping$mon_id <- c("MET","MSE","GLY")
with_chemistry <- ram_atlas_cohort_table(verified_with_chemistry,"P12345")
assert(identical(with_chemistry$mon_id[1:3],c("MET","MSE","GLY")),
  "Canonical CSV export must retain experimental modified residue chemistry.")
assert(nrow(cohort)==6L && length(unique(cohort$pdb_id))==2L,
       "Cohort should include only successfully verified entities.")
summary <- ram_atlas_cohort_summary(verified,"P12345")
assert(summary$verified_entities==2L && summary$observed_positions==4L &&
       summary$observed_entity_positions==6L,
       "Canonical cohort summary should not confuse entities with positions.")
support <- ram_atlas_cohort_position_support(verified,"P12345")
assert(identical(support$uniprot_resi,c(40L,41L,43L,44L)) &&
       identical(support$verified_entity_count,c(1L,2L,2L,1L)),
       "Observed canonical positions must be counted by verified entity.")
# Two PDB positions at the same author coordinate map to competing target
# positions: preserve both records but exclude that local residue from support.
bad <- verified
bad[["1ABC_1"]]$mapping <- rbind(
  bad[["1ABC_1"]]$mapping,
  transform(bad[["1ABC_1"]]$mapping[1L,,drop=FALSE],
            uniprot_resi=55L))
table <- ram_atlas_cohort_table(bad,"P12345")
assert(sum(table$mapping_status=="ambiguous")==2L,
       "Conflicting local SIFTS coordinates must remain ambiguous.")
assert(ram_atlas_cohort_summary(bad,"P12345")$ambiguous_local_residues==1L,
       "Ambiguous local residues must be counted once.")
assert(!55L %in% ram_atlas_cohort_position_support(bad,"P12345")$uniprot_resi,
       "Ambiguous mappings must not contribute to canonical support.")
wrong <- verified
wrong[["1ABC_1"]]$mapping$uniprot_accession[[1L]] <- "Q99999"
assert(inherits(try(ram_atlas_cohort_table(wrong,"P12345"),silent=TRUE),
                "try-error"),"Mixed UniProt cohorts must be rejected.")
message("Verified Atlas cohort coverage tests passed.")
