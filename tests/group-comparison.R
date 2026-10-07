# Run from repository root: Rscript tests/group-comparison.R
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


sparse_a <- prepared_a$models
sparse_a[[2L]] <- sparse_a[[2L]][sparse_a[[2L]]$resi!=1L,,drop=FALSE]
sparse_cmp <- ram_group_conformation_compare(
  reference,sparse_a,prepared_b$models,"apo","holo")
sparse_r1 <- sparse_cmp[sparse_cmp$resi==1L,,drop=FALSE]
assert(isTRUE(all.equal(sparse_r1$a_coverage,0.5)) &&
       identical(as.character(sparse_r1$evidence_profile),"Sparse coverage") &&
       !isTRUE(sparse_r1$high_support_shift),
       "Sparse group coverage must be explicit and must not be promoted as high support.")

bad <- with_decoy(make_chain("Y",c(0,0,0),c(0,0,0),"Favored"))
bad$resn[bad$chain=="Y"] <- c("GLY","GLY","GLY")
assert(inherits(try(ram_best_chain_match(reference,bad,
  min_identity=.8,min_reference_coverage=.8),silent=TRUE),"try-error"),
  "Unrelated candidate chains must fail explicit alignment criteria.")

message("Group conformation comparison tests passed.")
