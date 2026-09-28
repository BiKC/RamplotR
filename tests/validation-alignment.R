# Run from repository root: Rscript tests/validation-alignment.R
source(file.path("benchmarks", "validation_helpers.R"))
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
keys <- c("1.A.THR", "27.A.ALA", "27.A.ALA.B", "46.A.ASN")
reference <- c("  1.A.THR", "   27.A.ALA", "27.A.ALA.B", "46.A.ASN")
assert(identical(ram_match_reference_keys(keys, reference), seq_along(keys)),
       "Padded Bio3D keys and insertion codes must match unambiguously")
assert(is.na(ram_match_reference_keys("27.A.ALA.C", reference)),
       "A different insertion code must never be silently stripped")
ambiguous <- try(ram_match_reference_keys("27.A.ALA",
                        c("27.A.ALA", " 27.A.ALA")), silent = TRUE)
assert(inherits(ambiguous, "try-error"),
       "Ambiguous Bio3D identifiers must fail instead of silently matching")
message("Bio3D identifier alignment tests passed")
