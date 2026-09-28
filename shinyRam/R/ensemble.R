# Ensemble diagnostics. All comparisons use residue identifiers and circular
# angular statistics; an average of +179° and -179° is 180°, not zero.
# Does not infer prediction confidence or a single "best" model.

ram_ensemble_circular <- function(values) {
  z <- as.numeric(values[is.finite(values)])
  if(!length(z)) return(c(mean=NA_real_,sd=NA_real_))
  if(length(z)==1L) return(c(mean=z[[1L]],sd=NA_real_))
  radians <- z*pi/180
  x <- mean(cos(radians)); y <- mean(sin(radians))
  length_mean <- sqrt(x*x+y*y)
  mean_angle <- if(length_mean < 1e-8) NA_real_ else
    ((atan2(y,x)*180/pi+180) %% 360) - 180
  sd <- if(length_mean < 1e-8) NA_real_ else
    sqrt(-2*log(min(1,length_mean)))*180/pi
  c(mean=mean_angle,sd=sd)
}

ram_ensemble_key <- function(data) {
  paste(as.character(data$chain),data$resi,
        as.character(data$insertion_code),
        toupper(as.character(data$resn)),sep="\r")
}

ram_ensemble_summary <- function(models) {
  if(!is.list(models) || !length(models))
    stop("Provide at least one model's residue table.")
  required <- c("chain","resi","insertion_code","resn","phi","psi","region")
  for(tbl in models) {
    if(!is.data.frame(tbl) || !all(required %in% names(tbl)))
      stop("Every model needs identifiers, phi/psi angles and classifications.")
    if(anyDuplicated(ram_ensemble_key(tbl)))
      stop("Duplicate residue identifiers in an ensemble model.")
  }
  keys <- lapply(models,ram_ensemble_key)
  all_keys <- unique(unlist(keys,use.names=FALSE))
  n <- length(all_keys)
  nmodels <- length(models)
  if(!n) return(data.frame(chain=character(),resi=integer(),
      insertion_code=character(),resn=character(),models_present=integer(),
      phi_models=integer(),psi_models=integer(),
      phi_mean=numeric(),phi_sd=numeric(),psi_mean=numeric(),psi_sd=numeric(),
      classified_models=integer(),class_consistency=numeric(),
      region_mode=character(),changes_class=logical(),
      stringsAsFactors=FALSE))
  get <- function(col, type="numeric") {
    result <- matrix(if(type=="character") NA_character_ else NA_real_,
                     nrow=n,ncol=nmodels)
    for(i in seq_along(models)) {
      idx <- match(all_keys,keys[[i]])
      good <- which(!is.na(idx))
      if(length(good)) result[good,i] <- models[[i]][[col]][idx[good]]
    }
    result
  }
  phi <- get("phi");psi <- get("psi");region <- get("region","character")
  name_rows <- data.frame(chain=character(),resi=integer(),
    insertion_code=character(),resn=character(),stringsAsFactors=FALSE)
  # Use the first occurrence of each ID across all models, not row number.
  ids <- do.call(rbind,lapply(models,function(x)
    x[,required[1:4],drop=FALSE]))
  first <- !duplicated(unlist(keys,use.names=FALSE))
  info <- ids[first,,drop=FALSE]
  info <- info[match(all_keys,unlist(keys,use.names=FALSE)[first]),,drop=FALSE]
  calc <- function(matrix_values) {
    result <- t(vapply(seq_len(n),function(i)
      ram_ensemble_circular(matrix_values[i,]),numeric(2)))
    colnames(result) <- c("mean","sd")
    result
  }
  ph <- calc(phi);ps <- calc(psi)
  region_mode <- vapply(seq_len(n),function(i) {
    values <- region[i,]; values <- values[!is.na(values)]
    if(!length(values)) NA_character_ else
      names(sort(table(values),decreasing=TRUE))[[1L]]
  },character(1))
  agreement <- vapply(seq_len(n),function(i) {
    values <- region[i,];values <- values[!is.na(values)]
    if(!length(values)) NA_real_ else max(table(values))/length(values)
  },numeric(1))
  out <- data.frame(info,
    models_present=as.integer(rowSums(!is.na(phi)|!is.na(psi))),
    phi_models=as.integer(rowSums(is.finite(phi))),
    psi_models=as.integer(rowSums(is.finite(psi))),
    phi_mean=ph[,"mean"],phi_sd=ph[,"sd"],
    psi_mean=ps[,"mean"],psi_sd=ps[,"sd"],
    classified_models=as.integer(rowSums(!is.na(region))),
    class_consistency=as.numeric(agreement),region_mode=region_mode,
    changes_class=!is.na(agreement) & agreement<1,
    stringsAsFactors=FALSE,check.names=FALSE)
  out[order(-as.integer(out$changes_class),
            -pmax(replace(out$phi_sd,is.na(out$phi_sd),0),
                  replace(out$psi_sd,is.na(out$psi_sd),0)),
            out$chain,out$resi,out$insertion_code),,drop=FALSE]
}

ram_ensemble_analyze <- function(pdb,classifier,max_models=30L) {
  if(!exists("ram_model_count",mode="function") ||
     !exists("ram_model_at",mode="function") ||
     !exists("ram_extract_torsions",mode="function"))
    stop("Load structure IO and backbone functions first.")
  total <- ram_model_count(pdb)
  count <- min(total,as.integer(max_models))
  if(is.na(count) || count < 1L || count>100L)
    stop("Analyze between 1 and 100 structural models.")
  result <- lapply(seq_len(count),function(i) {
    single <- ram_model_at(pdb,i)
    torsions <- ram_extract_torsions(single)
    classifier(torsions)
  })
  list(summary=ram_ensemble_summary(result),analyzed_models=count,
       available_models=total,limited=count<total)
}
