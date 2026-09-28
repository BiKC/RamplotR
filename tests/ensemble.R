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
only_one <- other[other$resi!=1L,,drop=FALSE]
partial <- ram_ensemble_summary(list(m,only_one))
p <- partial[partial$resi==1L,,drop=FALSE]
assert(p$models_present==1L && p$phi_models==1L &&
       is.na(p$phi_sd),"A residue missing in model 2 must not be fabricated.")
bad <- rbind(m,m[1,,drop=FALSE])
assert(inherits(try(ram_ensemble_summary(list(bad)),silent=TRUE),
               "try-error"),"Ambiguous IDs must be rejected.")
message("Circular ensemble geometry and residue-alignment tests passed.")
