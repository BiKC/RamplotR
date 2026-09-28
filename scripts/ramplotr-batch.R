#!/usr/bin/env Rscript
# From the repository root:
# Rscript scripts/ramplotr-batch.R --input structures/ --output results/ \
#   --reference original --model 1 --report --ensemble-models 20
help <- paste(
  "RamplotR batch analysis (local structures; no network requests).",
  "Usage: Rscript scripts/ramplotr-batch.R --input FILE_OR_DIRECTORY --output DIR",
  "       [--reference original|alphafold|alphafold_filtered|astral2.08|custom_high_resolution]",
  "       [--background General|GLY|PRO|preProline] [--mode residue|legacy]",
  "       [--model N] [--ensemble-models N] [--max-files N]",
  "       [--prediction-source experimental|alphafold2|alphafold3|esmfold|other_prediction]",
  "       [--validation-xml PATH] [--report] [--no-json] [--overwrite]",
  "",
  "CSV output includes native diagnostic omega/chi1 angles. Rotamer/clash",
  "annotations are exclusively from the optional independent wwPDB XML.",
  "For ensemble analysis, the first N consistent atom-record models are analysed.",
  sep="\n"
)
required <- c("shinyRam/R/io.R","shinyRam/R/batch.R")
if(!all(file.exists(required)))
  stop("Run this script from the RamplotR repository root.",call.=FALSE)
for(filename in c("io","backbone","ramachandran","geometry","experimental",
                  "ensemble","predictions","reports","batch"))
  source(file.path("shinyRam","R",paste0(filename,".R")))
args <- commandArgs(trailingOnly=TRUE)
if("--help" %in% args) {cat(help,"\n");quit(status=0L)}
options <- tryCatch(ram_batch_options(args),error=function(e) {
  cat("Error:",conditionMessage(e),"\n\n",help,"\n",file=stderr())
  quit(status=2L)
})
if(options$json && !requireNamespace("jsonlite",quietly=TRUE))
  stop("Install jsonlite for machine-readable JSON, or supply --no-json.")
results <- ram_batch_run(options,repo_root=".")
print(results,row.names=FALSE)
quit(status=if(any(results$status=="error")) 1L else 0L)
