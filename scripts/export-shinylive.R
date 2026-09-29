# Run from repository root:
# Rscript scripts/export-shinylive.R bikc.be https://bikc.be/RamplotR/reference-data
# Normal Shiny and the original reference files are unchanged.
ramplotr_export_shinylive <- function(args) {
  if (length(args) != 2L || !grepl("^https?://", args[[2L]])) {
    stop("Usage: Rscript scripts/export-shinylive.R OUTPUT_DIR ABSOLUTE_REFERENCE_URL")
  }
  if (!requireNamespace("shinylive", quietly = TRUE)) {
    stop("Install the shinylive R package")
  }
  target <- normalizePath(args[[1L]], mustWork = FALSE)
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  base_url <- sub("/+$", "", args[[2L]])
  source_app <- normalizePath("shinyRam", mustWork = TRUE)
  datasets <- c("original", "alphafold", "alphafold_filtered",
                "astral2.08", "custom_high_resolution")
  stage <- tempfile("ramplotr-thin-")
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  thin <- file.path(stage, "shinyRam")
  dir.create(thin)
  dir.create(file.path(thin, "static"))

  # Keep scripts, assets and tiny reference indexes in app.json, not the RDS
  # grids. Ordinary Shiny uses the original complete shinyRam directory.
  for (entry in c("app.R", "R", "www")) {
    if (!file.copy(file.path(source_app, entry), thin, recursive = TRUE)) {
      stop("Failed to copy ", entry)
    }
  }
  staging <- tempfile("ramplotr-public-data-")
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  public <- file.path(staging, "reference-data")
  dir.create(public)

  for (dataset in datasets) {
    source_dir <- file.path(source_app, "static", dataset)
    names <- list.files(source_dir)
    names <- names[!file.info(file.path(source_dir, names))$isdir]
    if (!length(names)) stop("Missing reference dataset: ", dataset)
    from <- file.path(source_dir, names)
    manifest <- data.frame(
      file = names,
      md5 = unname(tools::md5sum(from)),
      stringsAsFactors = FALSE
    )
    dir.create(file.path(thin, "static", dataset))
    utils::write.table(
      manifest, file.path(thin, "static", dataset, "reference-index.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
    dest_dir <- file.path(public, dataset)
    dir.create(dest_dir)
    if (!all(file.copy(from, dest_dir))) {
      stop("Failed to copy reference files: ", dataset)
    }
    if (!identical(
      unname(tools::md5sum(file.path(dest_dir, names))), manifest$md5
    )) stop("Published references differ: ", dataset)
  }

  # This option is only added to the temporary app for the browser.
  app <- file.path(thin, "app.R")
  lines <- readLines(app, warn = FALSE)
  lines <- append(
    lines,
    paste0("options(ramplotr.reference_base_url = ", deparse(base_url), ")"),
    after = 0L
  )
  writeLines(lines, app)

  # Replacing only RamplotR retains shared shinylive assets and any other
  # static apps previously published beneath the same destination.
  destination <- file.path(target, "RamplotR")
  unlink(destination, recursive = TRUE)
  shinylive::export(
    thin, target, subdir = "RamplotR",
    template_params = list(
      title = "RamplotR | Protein structure analysis",
      include_in_head = paste0(
        '<link rel="icon" type="image/svg+xml" href="./favicon.svg">'
      )
    )
  )

  # The outer Shinylive page and Shiny's iframe both get the same favicon.
  if (!file.copy(
    file.path(source_app, "www", "favicon.svg"),
    file.path(destination, "favicon.svg"),
    overwrite = TRUE
  )) stop("Could not publish favicon")
  external <- file.path(destination, "reference-data")
  if (!file.rename(public, external)) {
    dir.create(external, recursive = TRUE, showWarnings = FALSE)
    if (!all(file.copy(list.files(public, full.names = TRUE),
                       external, recursive = TRUE))) {
      stop("Failed to publish reference-data")
    }
  }
  # Optional Apache settings for one.com. Scope rules to RamplotR and the
  # shared Shinylive asset folder. Never touch the website root .htaccess.
  # An existing shared .htaccess is left intact for manual merging.
  config <- file.path("config", "onecom")
  app_htaccess <- file.path(destination, ".htaccess")
  if (!file.copy(file.path(config, "ramplotr.htaccess"), app_htaccess,
                 overwrite = TRUE)) stop("Failed to include app .htaccess")
  shared_htaccess <- file.path(target, "shinylive", ".htaccess")
  if (!file.exists(shared_htaccess)) {
    if (!file.copy(file.path(config, "shinylive.htaccess"),
                   shared_htaccess)) {
      warning("Could not include the shared Shinylive .htaccess")
    }
  } else {
    message("Preserved existing shinylive/.htaccess; review it for cache rules.")
  }
  message("Thin Shinylive export complete: ", destination)
  message("Serve reference-data from: ", base_url)
  invisible(destination)
}

ramplotr_export_shinylive(commandArgs(trailingOnly = TRUE))
