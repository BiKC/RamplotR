# Reference access for normal and thin Shinylive deployments.
# A thin export keeps an index and fetches each original RDS file only when used.

ram_reference_index <- function(directory) {
  index_file <- file.path(directory, "reference-index.tsv")
  if (!file.exists(index_file)) return(NULL)
  index <- utils::read.delim(index_file, colClasses = "character",
                             stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(names(index), c("file", "md5")) || !nrow(index) ||
      anyDuplicated(index$file) ||
      any(!grepl("^[A-Za-z0-9_.-]+$", index$file)) ||
      any(!grepl("^[a-fA-F0-9]{32}$", index$md5))) {
    stop("Invalid reference index", call. = FALSE)
  }
  index
}

ram_reference_choices <- function(directory) {
  index <- ram_reference_index(directory)
  if (is.null(index)) list.files(directory) else index$file
}

ram_ensure_reference <- function(path) {
  if (file.exists(path)) return(invisible(path))
  index <- ram_reference_index(dirname(path))
  name <- basename(path)
  position <- if (is.null(index)) NA_integer_ else match(name, index$file)
  if (is.na(position)) stop("Reference file not found: ", path, call. = FALSE)
  base <- getOption("ramplotr.reference_base_url", "")
  if (length(base) != 1L || is.na(base) || !grepl("^https?://", base)) {
    stop("No reference download URL configured", call. = FALSE)
  }
  url <- paste0(sub("/+$", "", base), "/",
                basename(dirname(path)), "/", name)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile("reference-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  result <- tryCatch(
    utils::download.file(url, tmp, mode = "wb", quiet = TRUE),
    error = function(e) stop("Reference download failed: ",
                            conditionMessage(e), call. = FALSE)
  )
  if (!identical(result, 0L) || !file.exists(tmp)) {
    stop("Reference download failed: ", url, call. = FALSE)
  }
  actual <- unname(tools::md5sum(tmp))
  if (!identical(tolower(actual), tolower(index$md5[[position]]))) {
    stop("Reference checksum mismatch: ", name,
         ". Re-export the app and its reference-data together.", call. = FALSE)
  }
  if (!file.rename(tmp, path) && !file.exists(path)) {
    stop("Could not cache reference file: ", path, call. = FALSE)
  }
  invisible(path)
}
