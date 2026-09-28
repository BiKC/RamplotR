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
message("Publication SVG, PNG and self-contained HTML report tests passed.")
