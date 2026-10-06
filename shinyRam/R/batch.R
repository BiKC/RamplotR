# Local, reproducible batch analysis. No external service is contacted by default.
# The CLI entry point is scripts/ramplotr-batch.R.
ram_batch_options <- function(args) {
  keys <- c("input","output","reference","background","mode","model",
            "max-files","prediction-source","validation-xml","ensemble-models")
  flags <- c("report","no-json","overwrite","help")
  out <- list(input=NULL,output=NULL,reference="original",
    background="General",mode="residue",model=1L,max_files=100L,
    prediction_source="experimental",validation_xml=NULL,
    ensemble_models=0L,report=FALSE,json=TRUE,overwrite=FALSE,help=FALSE)
  i <- 1L
  while(i<=length(args)) {
    flag <- sub("^--","",args[[i]])
    if (!startsWith(args[[i]],"--") || !flag %in% c(keys,flags))
      stop("Unknown batch option: ",args[[i]],call.=FALSE)
    if(flag %in% flags) {
      name <- switch(flag,"no-json"="json",flag)
      out[[name]] <- flag!="no-json"
    } else {
      if(i==length(args) || startsWith(args[[i+1L]],"--"))
        stop("Missing value for --",flag,call.=FALSE)
      value <- args[[i+1L]]
      name <- gsub("-","_",flag,fixed=TRUE)
      out[[name]] <- value
      i <- i+1L
    }
    i <- i+1L
  }
  if(out$help) return(out)
  if(is.null(out$input) || is.null(out$output))
    stop("Both --input and --output are required. See --help.")
  if(!out$reference %in% c("original","alphafold","alphafold_filtered",
                            "astral2.08","custom_high_resolution"))
    stop("Unknown reference distribution.")
  if(!out$background %in% c("General","GLY","PRO","preProline"))
    stop("Unknown plotting background.")
  if(!out$mode %in% c("residue","legacy")) stop("Invalid classification mode.")
  if(!out$prediction_source %in% c("experimental","alphafold_db",
                 "alphafold2","alphafold3","esmfold","other_prediction"))
    stop("Unknown prediction provenance.")
  for(field in c("model","max_files","ensemble_models")) {
    value <- suppressWarnings(as.integer(out[[field]]))
    lower <- if(field=="ensemble_models") 0L else 1L
    upper <- if(field=="max_files") 1000L else 100L
    if(length(value)!=1L || is.na(value) || value<lower || value>upper)
      stop("Invalid value for ",field)
    out[[field]] <- value
  }
  if(!is.null(out$validation_xml) && dir.exists(out$input))
    stop("--validation-xml requires one structure file to avoid ambiguous matching.")
  out
}

ram_batch_files <- function(input,max_files=100L) {
  if(dir.exists(input)) {
    files <- list.files(input,pattern="\\.(pdb|ent|cif|mmcif|mcif)$",
                        full.names=TRUE,ignore.case=TRUE,recursive=FALSE)
  } else if(file.exists(input)) {
    ram_detect_format(input)
    files <- input
  } else stop("Input path does not exist.")
  files <- sort(files)
  if(!length(files) || length(files)>max_files)
    stop("No supported structures, or more than ",max_files," input files.")
  normalizePath(files,mustWork=TRUE)
}

ram_batch_run <- function(options,repo_root=".") {
  if(!requireNamespace("bio3d",quietly=TRUE))
    stop("Install Bio3D before running the batch analysis.")
  files <- ram_batch_files(options$input,options$max_files)
  destination <- normalizePath(options$output,mustWork=FALSE)
  dir.create(destination,recursive=TRUE,showWarnings=FALSE)
  if(!dir.exists(destination)) stop("Unable to create output directory.")
  refs <- file.path(repo_root,"shinyRam","static",options$reference)
  reference_file <- file.path(refs,options$background)
  background <- ram_read_reference(reference_file)
  official <- if(!is.null(options$validation_xml))
    ram_external_validation_read(options$validation_xml) else NULL
  source_names <- make.unique(gsub("[^A-Za-z0-9._-]","_",
    basename(files)))
  summary <- vector("list",length(files))
  for(i in seq_along(files)) {
    input <- files[[i]]
    prefix <- source_names[[i]]
    summary[[i]] <- tryCatch({
      base <- file.path(destination,prefix)
      targets <- c(paste0(base,".residues.csv"),
                   if(options$json) paste0(base,".json"),
                   if(options$report) c(paste0(base,".svg"),
                                       paste0(base,".html")),
                   if(options$ensemble_models>0L) paste0(base,".ensemble.csv"))
      if(!options$overwrite && any(file.exists(targets)))
        stop("Output already exists; use --overwrite to replace it.")
      pdb <- ram_load_structure(path=input,original_name=basename(input))
      selected <- ram_model_at(pdb,options$model)
      backbone <- ram_extract_torsions(selected)
      classified <- ram_classify_torsions(backbone,refs,background,
        options$mode,threshold_fn=ram_density_thresholds)
      classified <- ram_rama8000_classify(
        classified, file.path(repo_root,"shinyRam","static","rama8000"))
      diagnostic <- ram_extra_geometry(selected,backbone)
      classified <- ram_join_geometry(classified,diagnostic)
      if(options$prediction_source!="experimental") {
        confidence <- ram_prediction_from_atoms(selected,backbone,
                                                  options$prediction_source)
        classified$plddt <- confidence$plddt
        classified$confidence_category <- confidence$confidence_category
      }
      if(!is.null(official)) classified <-
        ram_external_validation_join(classified,official,options$model)
      metadata <- ram_report_metadata(basename(input),options$reference,
        options$background,options$mode,options$model,reference_file)
      metadata$input_md5 <- unname(tools::md5sum(input))
      metadata$prediction_provenance <- options$prediction_source
      metadata$native_diagnostics <- "peptide omega and descriptive chi1"
      metadata$rama8000_validation <- "six-class cctbx/Phenix-compatible evaluation"
      metadata$independent_wwpdb <- if(is.null(official)) "none" else
        basename(options$validation_xml)
      if(!is.null(official)) metadata$independent_wwpdb_md5 <-
        unname(tools::md5sum(options$validation_xml))
      utils::write.csv(classified,paste0(base,".residues.csv"),
                       row.names=FALSE,na="")
      ens <- NULL
      if(options$ensemble_models>0L && ram_model_count(pdb)>1L) {
        ens <- ram_ensemble_analyze(pdb,classifier=function(torsions)
          ram_classify_torsions(torsions,refs,background,options$mode,
                                threshold_fn=ram_density_thresholds),
          max_models=options$ensemble_models)
        utils::write.csv(ens$summary,paste0(base,".ensemble.csv"),
                         row.names=FALSE,na="")
        metadata$ensemble_analyzed <- ens$analyzed_models
        metadata$ensemble_available <- ens$available_models
      }
      if(options$json) {
        if(!requireNamespace("jsonlite",quietly=TRUE))
          stop("Install jsonlite or use --no-json.")
        jsonlite::write_json(list(metadata=metadata,
          residues=classified,ensemble=if(is.null(ens)) NULL else ens$summary,
          independent_summary=ram_external_validation_summary(classified)),
          paste0(base,".json"),pretty=TRUE,auto_unbox=TRUE,na="null",
          dataframe="rows",digits=8)
      }
      if(options$report) {
        svg <- paste0(base,".svg")
        ram_save_figure(svg,classified,background,
          c("#f1faf8","#b7d6ce","#69ada4","#126c73"),
          chain_colors=c(),format="svg",title=basename(input))
        ram_save_html_report(paste0(base,".html"),classified,metadata,svg)
      }
      data.frame(input=basename(input),status="ok",
        atoms=nrow(selected$atom),residues=nrow(classified),
        classified=sum(!is.na(classified$region)),error="",stringsAsFactors=FALSE)
    },error=function(e) data.frame(input=basename(input),status="error",
        atoms=NA_integer_,residues=NA_integer_,classified=NA_integer_,
        error=conditionMessage(e),stringsAsFactors=FALSE))
  }
  tab <- do.call(rbind,summary)
  utils::write.csv(tab,file.path(destination,"batch-summary.csv"),
                   row.names=FALSE)
  if(any(tab$status=="error"))
    warning("Some structures failed; inspect batch-summary.csv.")
  tab
}
