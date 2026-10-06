source("shinyRam/R/ensemble.R")
assert <- function(x,msg) if(!isTRUE(x)) stop(msg,call.=FALSE)
mean <- ram_ensemble_circular(c(179,-179))
assert(abs(abs(mean[["mean"]])-180)<2 && mean[["sd"]]<2,
       "Circular angle statistics must not average across the -180/180 seam.")
assert(is.na(ram_ensemble_circular(c(NA_real_))[["mean"]]) &&
       is.na(ram_ensemble_circular(c(179))[["sd"]]),
       "Missing observations and one-model uncertainty must remain explicit.")
m <- data.frame(chain="A",resi=c(1L,2L),insertion_code="",
  resn=c("SER","PRO"),phi=c(179,-60),psi=c(-179,145),
  region=c("Favoured","Allowed"),stringsAsFactors=FALSE)
other <- m[2:1,,drop=FALSE]
other$phi[other$resi==1L] <- -179
other$psi[other$resi==1L] <- 179
other$region[other$resi==2L] <- "Not allowed"
out <- ram_ensemble_summary(list(m,other))
one <- out[out$resi==1L,,drop=FALSE]
two <- out[out$resi==2L,,drop=FALSE]
assert(one$phi_models==2L && one$psi_models==2L &&
       abs(abs(one$phi_mean)-180)<2 && one$phi_sd<2 &&
       !one$changes_class,"Ensemble must match rows by residue identity.")
assert(two$changes_class && two$class_consistency==0.5 &&
       two$classified_models==2L,
       "Classification changes across models must be explicit.")

m$rama8000_region <- c("Favored","Allowed")
m$plddt <- c(95,88)
other$rama8000_region <- ifelse(other$resi==1L,"Outlier","Allowed")
other$plddt <- ifelse(other$resi==1L,91,92)
extended <- ram_ensemble_summary(list(m,other))
e1 <- extended[extended$resi==1L,,drop=FALSE]
e2 <- extended[extended$resi==2L,,drop=FALSE]
assert(e1$rama8000_changes && e1$rama8000_consistency==0.5 &&
       e1$rama8000_models==2L,
       "Prediction ensembles must expose Rama8000 category disagreement.")
assert(!e2$rama8000_changes && e2$rama8000_consistency==1,
       "Stable Rama8000 categories must remain explicit.")
assert(e1$plddt_models==2L && abs(e1$plddt_mean-93)<1e-12 &&
       e1$plddt_min==91 && e1$plddt_max==95,
       "Prediction ensemble pLDDT summary changed.")
only_one <- other[other$resi!=1L,,drop=FALSE]
partial <- ram_ensemble_summary(list(m,only_one))
p <- partial[partial$resi==1L,,drop=FALSE]
assert(p$models_present==1L && p$phi_models==1L &&
       is.na(p$phi_sd),"A residue missing in model 2 must not be fabricated.")
bad <- rbind(m,m[1,,drop=FALSE])
assert(inherits(try(ram_ensemble_summary(list(bad)),silent=TRUE),
               "try-error"),"Ambiguous IDs must be rejected.")
# Prediction-ensemble orchestration is tested with small stubs so this unit test
# remains independent of Bio3D and file I/O.
ram_model_at <- function(pdb,index=1L) pdb
ram_extract_torsions <- function(pdb) pdb$torsions
ram_prediction_key <- function(chain,resi,insertion_code="") {
  paste(chain,resi,insertion_code,sep="\r")
}
ram_plddt_category <- function(value) {
  ifelse(value>=90,"Very high",ifelse(value>=70,"Confident",
    ifelse(value>=50,"Low","Very low")))
}
ram_prediction_from_atoms <- function(pdb,torsions,source) {
  data.frame(chain=torsions$chain,resi=torsions$resi,
    insertion_code=torsions$insertion_code,
    plddt=pdb$plddt,stringsAsFactors=FALSE)
}
make_prediction <- function(phi_shift=0,plddt=c(92,81)) list(
  torsions=data.frame(chain="A",resi=1:2,insertion_code="",
    resn=c("ALA","SER"),phi=c(-60,-80)+phi_shift,
    psi=c(-45,150),stringsAsFactors=FALSE),
  plddt=plddt
)
prediction_classifier <- function(torsions) {
  torsions$region <- c("Favoured","Allowed")
  torsions$rama8000_region <- if(torsions$phi[[1L]] < -55)
    c("Favored","Allowed") else c("Allowed","Allowed")
  torsions
}
prediction <- ram_prediction_ensemble_analyze(
  list(make_prediction(0,c(94,80)),make_prediction(12,c(88,84))),
  classifier=prediction_classifier,source="alphafold2",
  labels=c("seed-1","seed-2")
)
assert(prediction$analyzed_models==2L && prediction$available_models==2L &&
       prediction$common_residues==2L,
       "Prediction ensemble model/residue coverage changed.")
assert(identical(prediction$labels,c("seed-1","seed-2")) &&
       identical(as.character(prediction$model_summary$model),
                 c("seed-1","seed-2")),
       "Prediction ensemble model labels must be preserved.")
p1 <- prediction$summary[prediction$summary$resi==1L,,drop=FALSE]
assert(p1$plddt_models==2L && abs(p1$plddt_mean-91)<1e-12 &&
       p1$rama8000_changes && p1$rama8000_consistency==0.5,
       "Prediction ensemble must combine pLDDT and Rama8000 disagreement.")
assert(abs(prediction$model_summary$plddt_mean[[1L]]-87)<1e-12 &&
       prediction$model_summary$residues[[1L]]==2L,
       "Per-model prediction provenance summary changed.")
assert(inherits(try(ram_prediction_ensemble_analyze(
  list(make_prediction(),make_prediction()),prediction_classifier,
  source="alphafold3"),silent=TRUE),"try-error"),
  "AF3 must not be silently interpreted as a B-factor prediction ensemble.")
assert(inherits(try(ram_prediction_ensemble_analyze(
  list(make_prediction()),prediction_classifier,
  source="esmfold"),silent=TRUE),"try-error"),
  "Prediction ensemble helper must require at least two models.")
assert(inherits(try(ram_prediction_ensemble_analyze(
  list(make_prediction(),make_prediction()),prediction_classifier,
  source="esmfold",max_models=31L),silent=TRUE),"try-error"),
  "Prediction ensemble helper must enforce the documented 30-model maximum.")

message("Circular ensemble geometry and residue-alignment tests passed.")
