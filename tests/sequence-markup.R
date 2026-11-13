# The vectorized renderer must retain residue controls and escape file-derived text.
source('shinyRam/R/inspection.R')
stopifnot(requireNamespace('htmltools',quietly=TRUE))
fixture <- data.frame(chain=c('A','A','B"<'),resi=c(10L,10L,20L),
  insertion_code=c('','A',''),resn=c('ALA','GLY','UNK'),
  region=c('Favoured',NA,'Not allowed'),plddt=c(95,NA,40),
  rama8000_region=c('Favored',NA,'Outlier'),rama8000_group='General')
chain <- ram_sequence_data(fixture)
markup <- as.character(ram_sequence_markup(chain,c(TRUE,FALSE,TRUE)))
stopifnot(
  length(gregexpr('<button ',markup,fixed=TRUE)[[1]])==3L,
  grepl('data-insertion="A"',markup,fixed=TRUE),
  grepl('disabled="disabled"',markup,fixed=TRUE),
  grepl('data-chain="B&quot;&lt;"',markup,fixed=TRUE),
  !grepl('data-chain="B"<',markup,fixed=TRUE),
  grepl('ram-seq-standard-missing',markup,fixed=TRUE),
  grepl('Missing angles',markup,fixed=TRUE),
  grepl('--ram-plddt-color:#126e74',markup,fixed=TRUE),
  grepl('data-plddt="40.0"',markup,fixed=TRUE),
  grepl('ram-seq-aa">X</span>',markup,fixed=TRUE),
  identical(as.character(ram_sequence_markup(chain[FALSE,],logical())),''))
message('Vectorized sequence controls, missing angles/confidence and HTML escaping passed.')
