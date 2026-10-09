source(file.path("shinyRam","R","canonical.R"))
source(file.path("shinyRam","R","atlas-geometry.R"))
source(file.path("shinyRam","R","atlas-construct.R"))
assert <- function(ok,message) if(!isTRUE(ok)) stop(message,call.=FALSE)

positions <- 1:40
fixture <- function(shift=100L,mon=rep("ALA",40L),observed=rep(TRUE,40L)) {
  map <- data.frame(chain="A",resi=positions+shift,insertion_code="",
    uniprot_accession="P12345",uniprot_resi=positions,entity_id=1L,
    struct_asym_id="X",label_seq_id=positions,observed=observed,
    mon_id=mon,canonical_source="PDBe updated mmCIF exact SIFTS",
    stringsAsFactors=FALSE)
  points <- lapply(positions,function(i) list(
    struct_asym_id="X",label_seq_id=i,
    x=i,y=sin(i),z=cos(i)))
  list(state="mapped",mapping=map,ca_points=points)
}
verified <- list("1ABC_1"=fixture(), "2XYZ_1"=fixture(900L))
ids <- names(verified)
audit <- ram_atlas_construct_audit(verified,"P12345",ids)
assert(!audit$differs && !audit$incomplete &&
       audit$pairs$common_observed==40L &&
       audit$pairs$chemistry_checked==40L &&
       audit$pairs$chemistry_differences==0L,
  "Identical mapped chemistry should pass construct audit despite shifted author numbering.")
assert(audit$pairs$unmatched_a==0L && audit$pairs$unmatched_b==0L,
  "Same canonical positions should have full observed overlap.")

variant <- verified
variant[["2XYZ_1"]]$mapping$mon_id[[12L]] <- "MSE"
review <- ram_atlas_construct_audit(variant,"P12345",ids)
assert(review$differs && review$pairs$chemistry_differences==1L &&
       grepl("12:ALA/MSE",review$pairs$difference_examples,fixed=TRUE),
  "Different experimental monomer identity must be shown with exact canonical position.")
assert(review$pairs$chemistry_checked==40L,
  "Modified polymer monomers are evidence, not automatically missing.")

unknown <- variant
unknown[["2XYZ_1"]]$mapping$mon_id[[12L]] <- NA_character_
unclear <- ram_atlas_construct_audit(unknown,"P12345",ids)
assert(!unclear$differs && unclear$incomplete &&
       unclear$pairs$chemistry_checked==39L &&
       unclear$pairs$chemistry_unknown==1L,
  "Missing chemistry must not be interpreted as agreement or a mutation.")
legacy <- verified
legacy[["2XYZ_1"]]$mapping$mon_id <- NULL
legacy_review <- ram_atlas_construct_audit(legacy,"P12345",ids)
assert(legacy_review$incomplete &&
       legacy_review$pairs$chemistry_unknown==40L,
  "Historical mapping without monomer field remains explicitly unassessed.")

partial <- verified
partial[["2XYZ_1"]]$mapping$observed[1:8] <- FALSE
subset <- ram_atlas_construct_audit(partial,"P12345",ids)
assert(subset$incomplete && subset$pairs$common_observed==32L &&
       subset$pairs$unmatched_a==8L &&
       subset$pairs$unmatched_b==0L &&
       subset$pairs$chemistry_checked==32L,
  "Missing observed construct regions should be reported separately from chemistry differences.")

conflicted <- verified
conflicted[["2XYZ_1"]]$mapping <- rbind(
  conflicted[["2XYZ_1"]]$mapping,
  transform(conflicted[["2XYZ_1"]]$mapping[12L,,drop=FALSE],
    uniprot_resi=12L,mon_id="MSE"))
unresolved <- ram_atlas_construct_audit(conflicted,"P12345",ids)
assert(!unresolved$differs && unresolved$pairs$common_observed==39L &&
       unresolved$pairs$unmatched_a==1L && unresolved$incomplete,
  "Contradictory exact mapping must exclude the residue, not choose a monomer arbitrarily.")

isoform <- verified
isoform[["2XYZ_1"]]$mapping$uniprot_accession <- "P12345-2"
assert(inherits(try(ram_atlas_construct_audit(isoform,"P12345",ids),
  silent=TRUE),"try-error"),
  "Different UniProt isoforms cannot be silently equated.")

multi <- c(verified,list("3ABC_1"=fixture(200L)))
multi[["3ABC_1"]]$mapping$mon_id[[15L]] <- "SER"
all_pairs <- ram_atlas_construct_audit(multi,"P12345",names(multi))
assert(all_pairs$pair_count==3L &&
       sum(all_pairs$pairs$chemistry_differences)==2L,
  "Entire selected cohort should be checked pairwise, not just the first two entries.")

message("Experimental Atlas construct and chemistry audit tests passed.")
