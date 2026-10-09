# Run from repository root: Rscript tests/group-comparison.R
source(file.path("shinyRam","R","conformation.R"))
source(file.path("shinyRam","R","inspection.R"))
source(file.path("shinyRam","R","ensemble.R"))
source(file.path("shinyRam","R","group-comparison.R"))
assert <- function(x,msg) if(!isTRUE(x)) stop(msg,call.=FALSE)

make_chain <- function(chain, phi, psi, rama="Favored") {
  data.frame(
    chain=chain,resi=1:3,insertion_code="",
    resn=c("ALA","SER","THR"),
    phi=phi,psi=psi,
    region=c("Favoured","Favoured","Favoured"),
    rama8000_region=rep(rama,3),
    stringsAsFactors=FALSE
  )
}
reference <- make_chain("A",c(170,-60,-70),c(10,-40,140))

with_decoy <- function(target, decoy_chain="Z") {
  decoy <- data.frame(
    chain=decoy_chain,resi=1:3,insertion_code="",
    resn=c("GLY","GLY","GLY"),
    phi=c(-30,-30,-30),psi=c(30,30,30),
    region="Allowed",rama8000_region="Allowed",
    stringsAsFactors=FALSE
  )
  rbind(target,decoy)
}

a1 <- with_decoy(make_chain("X",c(170,-60,-70),c(10,-40,140),"Favored"))
a2 <- with_decoy(make_chain("X",c(-170,-62,-72),c(12,-42,142),"Favored"))
b1 <- with_decoy(make_chain("Q",c(140,-61,-69),c(10,-39,141),"Allowed"))
b2 <- with_decoy(make_chain("Q",c(150,-59,-71),c(11,-41,139),"Allowed"))

prepared_a <- ram_prepare_structure_group(reference,list(a1,a2),c("apo-1","apo-2"))
prepared_b <- ram_prepare_structure_group(reference,list(b1,b2),c("holo-1","holo-2"))
assert(identical(as.character(prepared_a$model_summary$chain),c("X","X")),
       "Automatic chain matching should choose the sequence-compatible chain.")
assert(identical(as.character(prepared_b$model_summary$chain),c("Q","Q")),
       "Automatic chain matching changed for group B.")
assert(all(prepared_a$model_summary$identity==1) &&
       all(prepared_b$model_summary$identity==1),
       "Exact synthetic homologues should align at 100% identity.")

cmp <- ram_group_conformation_compare(
  reference,prepared_a$models,prepared_b$models,"apo","holo")
r1 <- cmp[cmp$resi==1L,,drop=FALSE]
assert(nrow(r1)==1L,"Reference residue should appear once in group comparison.")
assert(abs(abs(r1$a_phi_mean)-180)<1e-8,
       "Circular mean of +170/-170 must remain at the ±180° boundary.")
assert(isTRUE(all.equal(r1$b_phi_mean,145,tolerance=1e-8)),
       "Group B circular mean changed unexpectedly.")
assert(isTRUE(all.equal(r1$delta_phi,-35,tolerance=1e-8)),
       "Between-group angular difference must wrap correctly.")
assert(r1$angular_displacement>=34.9 && r1$angular_displacement<=35.1,
       "Combined group shift should reflect the wrapped phi difference.")
assert(isTRUE(r1$consistent_shift),
       "Low-dispersion groups with a >=30° shift should be flagged for navigation.")
assert(isTRUE(r1$high_support_shift) &&
       isTRUE(all.equal(r1$a_coverage,1)) &&
       isTRUE(all.equal(r1$b_coverage,1)) &&
       identical(as.character(r1$evidence_profile),"Low-dispersion shift"),
       "Fully covered low-dispersion shifts should receive a high-support evidence profile.")
assert(isTRUE(r1$rama8000_mode_changed),
       "Different standard-validation modes should remain visible.")

assert(all(c("a_basin_mode","b_basin_mode","basin_mode_changed") %in% names(cmp)),
       "Group comparisons must retain modal backbone states for both groups.")

state_reference <- make_chain("A",c(-63,-63,-63),c(-43,-43,-43))
state_a1 <- make_chain("A",c(-63,-63,-63),c(-43,-43,-43))
state_a2 <- make_chain("A",c(-65,-61,-64),c(-41,-45,-42))
state_b1 <- make_chain("A",c(-135,-63,-63),c(135,-43,-43))
state_b2 <- make_chain("A",c(-133,-62,-65),c(137,-44,-42))
state_cmp <- ram_group_conformation_compare(
  state_reference,list(state_a1,state_a2),list(state_b1,state_b2),"apo","holo")
state_r1 <- state_cmp[state_cmp$resi==1L,,drop=FALSE]
assert(state_r1$a_basin_mode=="Alpha-R" && state_r1$b_basin_mode=="Beta" &&
       isTRUE(state_r1$basin_mode_changed),
       "Group comparison must flag a modal backbone-state transition.")


sparse_a <- prepared_a$models
sparse_a[[2L]] <- sparse_a[[2L]][sparse_a[[2L]]$resi!=1L,,drop=FALSE]
sparse_cmp <- ram_group_conformation_compare(
  reference,sparse_a,prepared_b$models,"apo","holo")
sparse_r1 <- sparse_cmp[sparse_cmp$resi==1L,,drop=FALSE]
assert(isTRUE(all.equal(sparse_r1$a_coverage,0.5)) &&
       identical(as.character(sparse_r1$evidence_profile),"Sparse coverage") &&
       !isTRUE(sparse_r1$high_support_shift),
       "Sparse group coverage must be explicit and must not be promoted as high support.")

# Separate phi and psi counts must never inflate complete-pair support.
complement_a1 <- make_chain("A",c(-60,-60,-60),c(-40,-40,-40))
complement_a2 <- complement_a1
complement_a3 <- complement_a1
complement_a4 <- complement_a1
complement_a1$psi[1] <- NA_real_
complement_a2$phi[1] <- NA_real_
complement_a3$psi[1] <- NA_real_
complement_b <- lapply(seq_len(4),function(i)
  make_chain("A",c(-120,-60,-60),c(100,-40,-40)))
paired_cmp <- ram_group_conformation_compare(reference,
  list(complement_a1,complement_a2,complement_a3,complement_a4),
  complement_b,"A","B")
paired_r1 <- paired_cmp[paired_cmp$resi==1L,,drop=FALSE]
assert(paired_r1$a_phi_models==2L && paired_r1$a_psi_models==3L &&
       paired_r1$a_paired_angle_models==1L &&
       isTRUE(all.equal(paired_r1$a_coverage,0.25)) &&
       !isTRUE(paired_r1$high_support_shift),
       "Phi/psi from different members must not count as paired support.")
assert(identical(as.character(paired_r1$evidence_profile),
  "Unreplicated structural difference"),
  "Single paired observations cannot establish within-group consistency.")
single_cmp <- ram_group_conformation_compare(reference,
  list(make_chain("A",c(-60,-60,-60),c(-40,-40,-40))),
  list(make_chain("A",c(-125,-60,-60),c(100,-40,-40))),"A","B")
assert(all(single_cmp$evidence_profile==
  "Unreplicated structural difference") &&
  !any(single_cmp$high_support_shift),
  "One-versus-one comparisons must be explicitly labelled unreplicated.")

bad <- with_decoy(make_chain("Y",c(0,0,0),c(0,0,0),"Favored"))
bad$resn[bad$chain=="Y"] <- c("GLY","GLY","GLY")
assert(inherits(try(ram_best_chain_match(reference,bad,
  min_identity=.8,min_reference_coverage=.8),silent=TRUE),"try-error"),
  "Unrelated candidate chains must fail explicit alignment criteria.")

message("Group conformation comparison tests passed.")
