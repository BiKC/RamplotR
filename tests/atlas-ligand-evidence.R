source(file.path("shinyRam","R","atlas-ligand-evidence.R"))
assert <- function(x,message) if(!isTRUE(x))stop(message,call.=FALSE)
ids <- c("1ABC_1","2XYZ_1","3ABC_1")
geometry <- list(selected=ids,assignment=data.frame(
  entity=ids,geometric_group=c(1L,2L,2L)))
measured <- list(status="measured",total_nonwater_sites=2L,
  nearby_sites=1L,sites=list(list(comp_id="ATP",asym_id="B",
    auth_seq_id="900",min_distance_A=3.5,nearest_uniprot_resi=104L,
    heavy_atoms=16L)),warning="")
verified <- list("1ABC_1"=list(ligand_contacts=measured),
  "2XYZ_1"=list(ligand_contacts=list(status="measured",
    total_nonwater_sites=0L,sites=list(),warning="No non-water HETATM")),
  "3ABC_1"=list(ligand_contacts=list(status="unavailable",
    warning="No atom coordinates")))
res <- ram_atlas_observed_ligand_context(verified,geometry)
assert(res$measured==2L && res$unavailable==1L &&
       res$with_proximity==1L && nrow(res$contacts)==1L,
  "Measured proximity, no recorded contact and unavailable are distinct.")
assert(identical(res$contacts$nearest_uniprot_resi,104L) &&
       abs(res$contacts$min_distance_A-3.5)<1e-8,
  "Contact must retain the exact mapped UniProt residue and distance.")
assert(is.na(res$entries$nearby_nonwater_sites[[3L]]) &&
       res$entries$nearby_nonwater_sites[[2L]]==0L,
  "Missing coordinates are unknown, not ligand-free.")
bad <- verified
bad[["1ABC_1"]]$ligand_contacts$sites[[1L]]$min_distance_A <- 800
broken <- ram_atlas_observed_ligand_context(bad,geometry)
assert(broken$unavailable==2L &&
       !nrow(broken$contacts),
  "Invalid distances are rejected rather than published as proximity.")
assert(grepl("not",res$method,fixed=TRUE) &&
       !grepl("classified as apo",res$method,fixed=TRUE),
  "Observed contacts must not claim functional ligand states.")

# Per-site distances include all exact mapped residues when available.
verified_full <- verified
verified_full[["1ABC_1"]]$ligand_contacts$sites[[1L]]$residue_contacts <-
  list(list(uniprot_resi=104L,min_distance_A=3.5),
       list(uniprot_resi=105L,min_distance_A=4.2))
mapped <- ram_atlas_observed_ligand_context(verified_full,geometry)
assert(nrow(mapped$residue_contacts)==2L &&
       all(mapped$residue_contacts$scope=="all-mapped-contacts") &&
       identical(mapped$residue_contacts$uniprot_resi,c(104L,105L)),
  "All exact mapped residues near a component must be kept.")
assert(nrow(res$residue_contacts)==1L &&
       res$residue_contacts$scope[[1L]]=="nearest-only",
  "Old cached evidence must be identified as nearest-only, not comprehensive.")
# Both present and unknown records should remain distinct after invalid data.
invalid_second <- verified_full
invalid_second[["1ABC_1"]]$ligand_contacts$sites[[1L]]$residue_contacts[[2L]]$min_distance_A <- -1
reject <- ram_atlas_observed_ligand_context(invalid_second,geometry)
assert(reject$unavailable==2L && nrow(reject$residue_contacts)==0L,
  "Invalid per-residue distances must invalidate the whole entry.")
message("Observed ligand proximity scientific tests passed.")
