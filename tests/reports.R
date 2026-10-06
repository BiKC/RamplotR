# Executed in browser-preview CI after Shiny / Bio3D / htmltools are installed.
source(file.path("shinyRam","R","ramachandran.R"))
source(file.path("shinyRam","R","reports.R"))
data <- data.frame(
 chain=c("A","A"), resi=c(1L,2L), insertion_code=c("",""),
 resn=c("ALA","GLY"), phi=c(-60,NA_real_),psi=c(-45,NA_real_),
 region=c("Favoured",NA_character_),density=c(30,NA_real_)
)
axis <- seq(-180,180,by=30)
x <- outer(axis,axis,function(a,b) pmax(0,100-abs(a+60)-abs(b+45)))
ref <- list(x=axis,y=axis,z=t(x)+1)
cols <- c("#FFF8ED","#D4ECE7","#7DB9B5","#126E74")
base <- tempfile("ram-report")
svg <- paste0(base,".svg")
png <- paste0(base,".png")
html <- paste0(base,".html")
on.exit(unlink(c(svg,png,html)),add=TRUE)
ram_save_figure(svg,data,ref,cols,c(A="#CE6A4D"),format="svg")
stopifnot(file.exists(svg),file.info(svg)$size>1000L,
          any(grepl("<svg",readLines(svg,warn=FALSE),fixed=TRUE)))
ram_save_figure(png,data,ref,cols,c(A="#CE6A4D"),format="png")
stopifnot(file.exists(png),file.info(png)$size>1000L)
metadata <- list(structure="synthetic",reference_set="test",model="1")
ram_save_html_report(html,data,metadata,svg)
contents <- paste(readLines(html,warn=FALSE),collapse="\n")
stopifnot(grepl("RamplotR",contents,fixed=TRUE),
          grepl("synthetic",contents,fixed=TRUE),
          grepl("Missing angles",contents,fixed=TRUE),
          grepl("<svg",contents,fixed=TRUE))
ensemble_html <- paste0(base,"-prediction-ensemble.html")
on.exit(unlink(ensemble_html),add=TRUE)
ensemble_result <- list(
  summary=data.frame(
    chain=c("A","A"),resi=c(10L,11L),insertion_code=c("",""),
    resn=c("ALA","SER"),models_present=c(3L,3L),
    phi_models=c(3L,3L),psi_models=c(3L,3L),
    phi_mean=c(-60,-82),phi_sd=c(4.5,28.2),
    psi_mean=c(-45,151),psi_sd=c(3.8,14.1),
    rama8000_models=c(3L,3L),
    rama8000_mode=c("Favored","Allowed"),
    rama8000_consistency=c(1,2/3),
    rama8000_changes=c(FALSE,TRUE),
    plddt_models=c(3L,3L),plddt_mean=c(93,78),
    plddt_sd=c(2.1,12.4),plddt_min=c(90,62),plddt_max=c(96,89),
    stringsAsFactors=FALSE
  ),
  model_summary=data.frame(
    model=c("seed","seed","seed-3"),residues=rep(2L,3),
    finite_phi_psi=rep(2L,3),rama8000_outliers=c(0L,0L,1L),
    plddt_mean=c(87,85,82),plddt_min=c(78,75,62),
    stringsAsFactors=FALSE
  ),
  provenance=data.frame(
    model=c("seed","seed","seed-3"),source="esmfold",
    input_role=c("uploaded","uploaded","uploaded"),
    structure_model=rep(1L,3),
    coordinate_md5=c("aaa","bbb","ccc"),stringsAsFactors=FALSE
  ),
  common_residues=2L,analyzed_models=3L,available_models=3L,
  source="esmfold",limited=FALSE
)
ram_save_prediction_ensemble_report(
  ensemble_html,ensemble_result,
  metadata=list(structure="synthetic prediction",reference_set="test")
)
ensemble_contents <- paste(readLines(ensemble_html,warn=FALSE),collapse="\n")
stopifnot(file.exists(ensemble_html),file.info(ensemble_html)$size>1500L,
          grepl("Prediction ensemble",ensemble_contents,fixed=TRUE),
          grepl("prediction uncertainty or heterogeneity",ensemble_contents,
                fixed=TRUE),
          grepl("coordinate_md5",ensemble_contents,fixed=TRUE),
          grepl("input_role",ensemble_contents,fixed=TRUE),
          grepl("aaa",ensemble_contents,fixed=TRUE),
          grepl("bbb",ensemble_contents,fixed=TRUE),
          grepl("Rama8000 disagreements",ensemble_contents,fixed=TRUE),
          grepl("28.2",ensemble_contents,fixed=TRUE))

message("Publication SVG, PNG and self-contained HTML report tests passed.")
