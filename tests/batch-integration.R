# End-to-end offline batch regression with known downloaded wwPDB structures.
# Called by ui-preview.yml after fetching 1CRN and 1D3Z.
for(filename in c("io","backbone","ramachandran","geometry","experimental",
                  "ensemble","predictions","reports","batch"))
  source(file.path("shinyRam","R",paste0(filename,".R")))
assert <- function(x,msg) if(!isTRUE(x)) stop(msg,call.=FALSE)
root <- "benchmarks/output/ui-preview"
one <- file.path(root,"1CRN.pdb")
nmr <- file.path(root,"1D3Z.pdb")
if(!all(file.exists(c(one,nmr)))) stop("Download 1CRN and 1D3Z first.")
output <- tempfile("ram-phase-c-batch-")
on.exit <- NULL
a <- ram_batch_options(c("--input",one,"--output",output,"--report"))
run <- ram_batch_run(a)
prefix <- file.path(output,"1CRN.pdb")
assert(nrow(run)==1L && run$status[[1]]=="ok" &&
       run$residues[[1]]>=40 &&
       file.exists(paste0(prefix,".residues.csv")) &&
       file.exists(paste0(prefix,".json")) &&
       file.exists(paste0(prefix,".html")) &&
       file.exists(paste0(prefix,".svg")),
       paste("Batch figure/report/JSON export failed:",run$error[[1]]))
csv <- utils::read.csv(paste0(prefix,".residues.csv"))
assert(all(c("omega","omega_status","chi1","phi","psi","region") %in%
           names(csv)) &&
       all(csv$omega_status %in% c("Cis","Trans","Twisted","Missing")),
       "Batch must export new geometry without changing original columns.")
doc <- paste(readLines(paste0(prefix,".html"),warn=FALSE),collapse="\n")
assert(grepl("Additional descriptive geometry",doc,fixed=TRUE) &&
       grepl("Analysis provenance",doc,fixed=TRUE),
       "Batch HTML report must contain geometry and provenance.")
if(!requireNamespace("jsonlite",quietly=TRUE))
  stop("jsonlite is required by the batch integration workflow.")
json <- jsonlite::fromJSON(paste0(prefix,".json"))
assert(json$metadata$reference_set=="original" &&
       !is.null(json$metadata$input_md5),
       "Machine-readable output must preserve scientific provenance.")
b <- ram_batch_options(c("--input",nmr,"--output",output,
  "--ensemble-models","2","--no-json"))
run2 <- ram_batch_run(b)
assert(nrow(run2)==1L && run2$status[[1]]=="ok" &&
       file.exists(file.path(output,"1D3Z.pdb.ensemble.csv")),
       paste("Multi-model ensemble batch failed:",run2$error[[1]]))
ensemble <- utils::read.csv(file.path(output,"1D3Z.pdb.ensemble.csv"))
assert(nrow(ensemble)>=50 &&
       all(c("phi_mean","phi_sd","changes_class","models_present") %in%
         names(ensemble)),
       "NMR ensemble CSV is missing circular statistics.")
repeat_error <- ram_batch_run(a)
assert(repeat_error$status[[1]]=="error" &&
       grepl("already exists",repeat_error$error[[1]]),
       "Batch runs must not silently overwrite earlier published results.")
unlink(output,recursive=TRUE)
message("Real-structure batch figure, JSON, HTML and NMR-ensemble tests passed.")
