# Run with Rscript tests/reference-loader.R from the repository root.
source("shinyRam/R/reference-loader.R")
source("shinyRam/R/ramachandran.R")
root <- tempfile("ramplotr-reference-test-")
dir.create(root)
original <- file.path(root, "original")
dir.create(original)
grid <- list(x = -1:1, y = -1:1, z = matrix(1:9, nrow = 3))
saveRDS(grid, file.path(original, "General"), compress = "gzip")
stopifnot(identical(ram_read_reference(file.path(original, "General")), grid))
stopifnot("General" %in% ram_reference_choices(original))
manifest <- data.frame(file = "General",
                       md5 = unname(tools::md5sum(file.path(original, "General"))))
write.table(manifest, file.path(original, "reference-index.tsv"),
            row.names = FALSE, quote = FALSE, sep = "\t")
stopifnot(identical(ram_reference_choices(original), "General"))
stopifnot(identical(ram_reference_index(original)$md5, manifest$md5))
stopifnot(identical(ram_read_reference(file.path(original, "General")), grid))
writeLines("bad data", file.path(original, "reference-index.tsv"))
stopifnot(inherits(try(ram_reference_index(original), silent = TRUE),
                   "try-error"))
unlink(root, recursive = TRUE)
cat("Local reference-loader and manifest checks passed.\n")
