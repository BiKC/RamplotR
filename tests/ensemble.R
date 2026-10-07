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
  "AF3 ensembles must reject models without one-to-one confidence sidecars.")

ram_prepare_prediction <- function(pdb,torsions,source,sidecar=NULL,
                                   summary_file=NULL,model_id="") {
  residues <- data.frame(
    chain=torsions$chain,resi=torsions$resi,
    insertion_code=torsions$insertion_code,
    plddt=pdb$plddt,stringsAsFactors=FALSE
  )
  list(
    source=source,residues=residues,pae=NULL,pae_rows=integer(),
    ptm=pdb$ptm,iptm=pdb$iptm,ranking_score=pdb$ranking_score,
    fraction_disordered=pdb$fraction_disordered,has_clash=pdb$has_clash,
    notes=character(),model_id=model_id,
    confidence_file=basename(sidecar),
    summary_file=if(is.null(summary_file)) "" else basename(summary_file)
  )
}
af3_a <- make_prediction(0,c(96,82))
af3_b <- make_prediction(8,c(90,78))
af3_a$ptm <- 0.82; af3_a$iptm <- 0.71; af3_a$ranking_score <- 0.76
af3_a$fraction_disordered <- 0.10; af3_a$has_clash <- FALSE
af3_b$ptm <- 0.79; af3_b$iptm <- 0.68; af3_b$ranking_score <- 0.70
af3_b$fraction_disordered <- 0.14; af3_b$has_clash <- TRUE
sidecars <- c(tempfile(fileext="_confidences.json"),
              tempfile(fileext="_confidences.json"))
summaries <- c(tempfile(fileext="_summary_confidences.json"),
               tempfile(fileext="_summary_confidences.json"))
file.create(c(sidecars,summaries))
af3_ensemble <- ram_prediction_ensemble_analyze(
  list(af3_a,af3_b),prediction_classifier,source="alphafold3",
  labels=c("seed-1_sample-0","seed-1_sample-1"),
  sidecars=sidecars,summary_files=summaries
)
assert(af3_ensemble$analyzed_models==2L &&
       abs(af3_ensemble$summary$plddt_mean[
         af3_ensemble$summary$resi==1L]-93)<1e-12,
       "AF3 ensemble must aggregate pLDDT from matched model confidence data.")
assert(identical(as.character(af3_ensemble$model_summary$model),
                 c("seed-1_sample-0","seed-1_sample-1")) &&
       isTRUE(all.equal(af3_ensemble$model_summary$ranking_score,c(0.76,0.70))) &&
       identical(af3_ensemble$model_summary$has_clash,c(FALSE,TRUE)),
       "AF3 ensemble must preserve per-sample ranking and clash provenance.")
unlink(c(sidecars,summaries))
paired <- ram_af3_pair_files(
  model_names=c(
    "job_seed-7_sample-0_model.cif",
    "job_seed-7_sample-1_model.cif"
  ),
  model_paths=c("/tmp/model0.cif","/tmp/model1.cif"),
  confidence_names=c(
    "job_seed-7_sample-1_confidences.json",
    "job_seed-7_sample-0_confidences.json"
  ),
  confidence_paths=c("/tmp/conf1.json","/tmp/conf0.json"),
  summary_names=c("job_seed-7_sample-0_summary_confidences.json"),
  summary_paths=c("/tmp/summary0.json")
)
assert(identical(paired$label,
  c("job_seed-7_sample-0","job_seed-7_sample-1")) &&
  identical(paired$confidence_path,c("/tmp/conf0.json","/tmp/conf1.json")) &&
  identical(paired$summary_path,c("/tmp/summary0.json","")),
  "AF3 files must pair by seed/sample filename stem, not upload order.")
assert(inherits(try(ram_af3_pair_files(
  c("job_seed-1_sample-0_model.cif","job_seed-1_sample-1_model.cif"),
  c("/m0","/m1"),
  c("job_seed-1_sample-0_confidences.json"),
  c("/c0")),silent=TRUE),"try-error"),
  "AF3 pairing must reject a model without its full confidence sidecar.")
assert(inherits(try(ram_af3_pair_files(
  c("model0.cif","model1.cif"),c("/m0","/m1"),
  c("conf0.json","conf1.json"),c("/c0","/c1")),silent=TRUE),"try-error"),
  "AF3 pairing must reject filenames that cannot establish sample identity.")

assert(inherits(try(ram_prediction_ensemble_analyze(
  list(make_prediction()),prediction_classifier,
  source="esmfold"),silent=TRUE),"try-error"),
  "Prediction ensemble helper must require at least two models.")
assert(inherits(try(ram_prediction_ensemble_analyze(
  list(make_prediction(),make_prediction()),prediction_classifier,
  source="esmfold",max_models=31L),silent=TRUE),"try-error"),
  "Prediction ensemble helper must enforce the documented 30-model maximum.")

message("Circular ensemble geometry and residue-alignment tests passed.")
