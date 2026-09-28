source(file.path("shinyRam","R","io.R"))
source(file.path("shinyRam","R","backbone.R"))
source(file.path("shinyRam","R","inspection.R"))
name <- "benchmarks/output/ui-preview/1D3Z.pdb"
if (!file.exists(name)) stop("Download the 1D3Z NMR fixture first.")
pdb <- ram_load_structure(path=name,original_name="1D3Z.pdb")
count <- ram_model_count(pdb)
stopifnot(count > 1L)
first <- ram_model_at(pdb,1L)
last <- ram_model_at(pdb,count)
stopifnot(nrow(first$atom)==nrow(last$atom))
a <- ram_extract_torsions(first)
b <- ram_extract_torsions(last)
stopifnot(identical(a$resn,b$resn),
          identical(a$chain,b$chain),
          identical(a$resi,b$resi),
          nrow(a)>0L)
differences <- abs(ram_angular_difference(a$phi,b$phi))
if (!any(is.finite(differences) & differences > 0.01))
  stop("NMR model switching did not change torsion angles.")
message(sprintf("Validated %d independently selectable NMR models.",count))
