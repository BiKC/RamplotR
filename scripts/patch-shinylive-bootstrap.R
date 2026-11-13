# Keep the shared vendor runtime untouched; replace only this app's entry point.
ram_patch_shinylive_bootstrap <- function(destination, source = file.path('scripts','shinylive-bootstrap.mjs')) {
  index <- file.path(destination,'index.html')
  html <- paste(readLines(index,warn=FALSE),collapse='\n')
  fingerprint <- substr(unname(tools::md5sum(source)),1L,12L)
  name <- paste0('ramplotr-bootstrap-',fingerprint,'.js')
  if(!file.copy(source,file.path(destination,name),overwrite=TRUE)) stop('Could not stage browser bootstrap')
  # These are generated entry-point scripts, not scripts inside the Shiny app.
  html <- gsub('(?s)<script\\b[^>]*>.*?</script>','',html,perl=TRUE)
  html <- sub('</head>',paste0('<script type="module" src="./',name,'"></script>\n</head>'),html,fixed=TRUE)
  if(grepl('load-shinylive-sw.js',html,fixed=TRUE)||!grepl(name,html,fixed=TRUE)) stop('Could not replace browser entry point')
  writeLines(html,index)
  invisible(name)
}
