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
  membership <- matrix(FALSE,nrow=n,ncol=nmodels)
  for(i in seq_along(models))
    membership[,i] <- !is.na(match(all_keys,keys[[i]]))

  optional_all <- function(field)
    all(vapply(models,function(tbl) field %in% names(tbl),logical(1)))

  if(!n) {
    empty <- data.frame(chain=character(),resi=integer(),
      insertion_code=character(),resn=character(),models_present=integer(),
      phi_models=integer(),psi_models=integer(),
      phi_mean=numeric(),phi_sd=numeric(),psi_mean=numeric(),psi_sd=numeric(),
      classified_models=integer(),class_consistency=numeric(),
      region_mode=character(),changes_class=logical(),
      stringsAsFactors=FALSE)
    if(optional_all("rama8000_region")) {
      empty$rama8000_models <- integer()
      empty$rama8000_consistency <- numeric()
      empty$rama8000_mode <- character()
      empty$rama8000_changes <- logical()
    }
    if(optional_all("plddt")) {
      empty$plddt_models <- integer()
      empty$plddt_mean <- numeric()
      empty$plddt_sd <- numeric()
      empty$plddt_min <- numeric()
      empty$plddt_max <- numeric()
    }
    return(empty)
  }

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
  phi <- get("phi"); psi <- get("psi"); region <- get("region","character")

  ids <- do.call(rbind,lapply(models,function(x)
    x[,required[1:4],drop=FALSE]))
  flat_keys <- unlist(keys,use.names=FALSE)
  first <- !duplicated(flat_keys)
  info <- ids[first,,drop=FALSE]
  info <- info[match(all_keys,flat_keys[first]),,drop=FALSE]

  calc <- function(matrix_values) {
    result <- t(vapply(seq_len(n),function(i)
      ram_ensemble_circular(matrix_values[i,]),numeric(2)))
    colnames(result) <- c("mean","sd")
    result
  }
  mode_and_consistency <- function(matrix_values) {
    mode <- vapply(seq_len(n),function(i) {
      values <- matrix_values[i,]; values <- values[!is.na(values)]
      if(!length(values)) NA_character_ else
        names(sort(table(values),decreasing=TRUE))[[1L]]
    },character(1))
    consistency <- vapply(seq_len(n),function(i) {
      values <- matrix_values[i,]; values <- values[!is.na(values)]
      if(!length(values)) NA_real_ else max(table(values))/length(values)
    },numeric(1))
    list(mode=mode,consistency=consistency,
         count=as.integer(rowSums(!is.na(matrix_values))))
  }

  ph <- calc(phi); ps <- calc(psi)
  native <- mode_and_consistency(region)
  out <- data.frame(info,
    models_present=as.integer(rowSums(membership)),
    phi_models=as.integer(rowSums(is.finite(phi))),
    psi_models=as.integer(rowSums(is.finite(psi))),
    phi_mean=ph[,"mean"],phi_sd=ph[,"sd"],
    psi_mean=ps[,"mean"],psi_sd=ps[,"sd"],
    classified_models=native$count,
    class_consistency=as.numeric(native$consistency),
    region_mode=native$mode,
    changes_class=!is.na(native$consistency) & native$consistency<1,
    stringsAsFactors=FALSE,check.names=FALSE)

  if(optional_all("rama8000_region")) {
    standard <- get("rama8000_region","character")
    stat <- mode_and_consistency(standard)
    out$rama8000_models <- stat$count
    out$rama8000_consistency <- stat$consistency
    out$rama8000_mode <- stat$mode
    out$rama8000_changes <- !is.na(stat$consistency) & stat$consistency<1
  }

  if(optional_all("plddt")) {
    confidence <- get("plddt")
    finite_count <- rowSums(is.finite(confidence))
    safe_stat <- function(fun) vapply(seq_len(n),function(i) {
      values <- confidence[i,is.finite(confidence[i,])]
      if(!length(values)) NA_real_ else fun(values)
    },numeric(1))
    out$plddt_models <- as.integer(finite_count)
    out$plddt_mean <- safe_stat(mean)
    out$plddt_sd <- vapply(seq_len(n),function(i) {
      values <- confidence[i,is.finite(confidence[i,])]
      if(length(values)<2L) NA_real_ else stats::sd(values)
    },numeric(1))
    out$plddt_min <- safe_stat(min)
    out$plddt_max <- safe_stat(max)
  }

  spread <- pmax(replace(out$phi_sd,is.na(out$phi_sd),0),
                 replace(out$psi_sd,is.na(out$psi_sd),0))
  standard_change <- if("rama8000_changes" %in% names(out))
    as.integer(out$rama8000_changes) else rep(0L,nrow(out))
  out[order(-standard_change,-as.integer(out$changes_class),-spread,
            out$chain,out$resi,out$insertion_code),,drop=FALSE]
}

ram_prediction_ensemble_analyze <- function(pdbs, classifier, source,
                                            labels = NULL,
                                            max_models = 30L) {
  if(!is.list(pdbs) || !length(pdbs))
    stop("Upload at least one predicted structure.")
  permitted <- c("alphafold2","esmfold","other_prediction")
  if(length(source)!=1L || !source %in% permitted)
    stop("Prediction ensembles currently support AF2/ColabFold, ESMFold or other pLDDT-in-B-factor models.")
  count <- min(length(pdbs),as.integer(max_models))
  if(!is.finite(count) || count<1L || count>100L)
    stop("Analyze between 1 and 100 prediction models.")
  if(is.null(labels)) labels <- paste("Model",seq_along(pdbs))
  labels <- as.character(labels)
  if(length(labels)!=length(pdbs) || any(!nzchar(labels)))
    stop("Every prediction model needs a label.")

  models <- lapply(seq_len(count),function(i) {
    pdb <- ram_model_at(pdbs[[i]],1L)
    torsions <- ram_extract_torsions(pdb)
    classified <- classifier(torsions)
    confidence <- ram_prediction_from_atoms(pdb,torsions,source)
    keys <- ram_prediction_key(classified$chain,classified$resi,
                               classified$insertion_code)
    confidence_keys <- ram_prediction_key(confidence$chain,confidence$resi,
                                          confidence$insertion_code)
    if(anyDuplicated(confidence_keys))
      stop("Prediction residue identifiers must be unique.")
    classified$plddt <- confidence$plddt[match(keys,confidence_keys)]
    classified$confidence_category <- ram_plddt_category(classified$plddt)
    classified$model_label <- labels[[i]]
    classified
  })
  summary <- ram_ensemble_summary(models)
  model_summary <- do.call(rbind,lapply(seq_len(count),function(i) {
    table <- models[[i]]
    finite_angles <- is.finite(table$phi) & is.finite(table$psi)
    data.frame(
      model=labels[[i]],
      residues=nrow(table),
      finite_phi_psi=sum(finite_angles),
      rama8000_outliers=if ("rama8000_region" %in% names(table))
        sum(table$rama8000_region=="Outlier",na.rm=TRUE) else NA_integer_,
      plddt_mean=if ("plddt" %in% names(table) && any(is.finite(table$plddt)))
        mean(table$plddt[is.finite(table$plddt)]) else NA_real_,
      plddt_min=if ("plddt" %in% names(table) && any(is.finite(table$plddt)))
        min(table$plddt[is.finite(table$plddt)]) else NA_real_,
      stringsAsFactors=FALSE
    )
  }))
  list(
    summary=summary,
    model_summary=model_summary,
    models=models,
    labels=labels[seq_len(count)],
    analyzed_models=count,
    available_models=length(pdbs),
    common_residues=sum(summary$models_present==count),
    source=source,
    limited=count<length(pdbs)
  )
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
