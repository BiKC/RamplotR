# UI contract checks use base R and do not need to launch a Shiny server.
app <- paste(readLines(file.path("shinyRam", "app.R"), warn = FALSE),
             collapse = "\n")
ui <- strsplit(app, "# Define server logic required to draw a histogram",
               fixed = TRUE)[[1L]][[1L]]
css <- paste(readLines(file.path("shinyRam", "www", "styles.css"),
                       warn = FALSE), collapse = "\n")
js <- paste(readLines(file.path("shinyRam", "www", "custom.js"),
                      warn = FALSE), collapse = "\n")
required <- c(
  "PDB", "inputSource", "structfile", "submit", "validationMode",
  "bgtype", "background", "AA", "chains", "chainColors",
  "ligands", "dna", "rna", "spinning", "rocking", "colorscheme",
  "bg1", "bg2", "bg3", "bg4", "plotly", "NGL",
  "comparePlot", "compareSwap", "NGLCompare",
  "regionselect", "regions", "summary",
  "atlasAccession", "atlasDiscover", "atlasStatus", "atlasResults",
  "openGuide", "groupGuide", "atlasGuide", "groupComparisonExports"
)
missing <- required[!vapply(required, function(id) {
  grepl(paste0('"', id, '"'), ui, fixed = TRUE)
}, logical(1))]
if (length(missing)) stop("Missing original UI ID(s): ",
                          paste(missing, collapse = ", "))
stopifnot(
  grepl('href = "styles.css"', ui, fixed = TRUE),
  grepl('href = "favicon.svg"', ui, fixed = TRUE),
  grepl('src = "plotly-loader.js"', ui, fixed = TRUE),
  !grepl('tags$script(src = "https://cdn.plot.ly', ui, fixed = TRUE),
  grepl('src = "custom.js"', ui, fixed = TRUE),
  grepl('src = "canonical-mapping.js"', ui, fixed = TRUE),
  grepl('src = "atlas-discovery.js"', ui, fixed = TRUE),
  grepl('src = "atlas-sifts-exact.js"', ui, fixed = TRUE),
  grepl('src = "experimental-search.js"', ui, fixed = TRUE),
  grepl('class = "ram-workspace"', ui, fixed = TRUE),
  grepl('class = "ram-charts"', ui, fixed = TRUE),
  grepl("tabsetPanel(", ui, fixed = TRUE),
  grepl('source(file.path("R", "guide.R")', ui, fixed = TRUE),
  grepl('value = "guide"', ui, fixed = TRUE),
  grepl('.ram-guide-steps', css, fixed = TRUE),
  grepl('.ram-group-compare-controls', css, fixed = TRUE),
  grepl("@media (max-width: 580px)", css, fixed = TRUE),
  grepl(".ram-plot-empty[hidden]", css, fixed = TRUE),
  grepl("Plotly.react(", js, fixed = TRUE),
  grepl('scaleanchor: "x"', js, fixed = TRUE),
  !grepl("Plotly.purge(", js, fixed = TRUE),
  !grepl('style = "height:800px;width:891px"', ui, fixed = TRUE)
)
message("UI contract and responsive-layout checks passed.")
