source("shinyRam/R/io.R")
source("shinyRam/R/batch.R")
assert <- function(x,msg) if(!isTRUE(x)) stop(msg,call.=FALSE)
opts <- ram_batch_options(c("--input","/tmp/one.pdb","--output","/tmp/result",
  "--reference","original","--report","--ensemble-models","2",
  "--prediction-source","esmfold"))
assert(opts$model==1L && opts$ensemble_models==2L && opts$report &&
       opts$json && opts$prediction_source=="esmfold",
       "Batch argument parsing lost scientific defaults.")
assert(!ram_batch_options(c("--input","input.pdb","--output","out",
       "--no-json"))$json,"JSON must be optional.")
reject <- function(args) inherits(try(ram_batch_options(args),silent=TRUE),
                                  "try-error")
assert(reject(c("--input","a","--output","b","--mode","invented")),
       "Unknown classification modes must be rejected.")
assert(reject(c("--input","a","--output","b","--max-files","50000")),
       "Resource limits must be enforced.")
assert(reject(c("--input","a","--output","b","--prediction-source","experimental",
                "--validation-xml","x","--bad-option")),
       "Unknown options must be rejected.")
temp <- tempfile("ram-batch-");dir.create(temp)
file.create(file.path(temp,c("B.pdb","a.cif","notes.txt")))
files <- ram_batch_files(temp)
assert(identical(basename(files),sort(c("B.pdb","a.cif"))),
       "Only supported local structure files may enter a batch.")
assert(inherits(try(ram_batch_files(temp,max_files=1),silent=TRUE),
                "try-error"),"Max-file guard must apply to directories.")
unlink(temp,recursive=TRUE)
message("Batch CLI argument and local file-discovery tests passed.")
