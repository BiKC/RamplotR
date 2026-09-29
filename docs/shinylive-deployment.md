# Thin Shinylive deployment

The normal app and batch CLI continue using all five local RDS reference
directories. A thin browser deployment keeps the same reference options,
but moves the original, unchanged RDS files out of the initial app.json.

## Build

From the repository root, with shinylive installed:

```r
system2("Rscript", c("scripts/export-shinylive.R", "bikc.be",
  "https://bikc.be/RamplotR/reference-data"))
```

On Windows in an R session you can also use:

```r
source("scripts/export-shinylive.R", echo = TRUE)
```

The sourced form needs `commandArgs(trailingOnly = TRUE)`, so the
Rscript invocation above is the supported approach. For a custom hostname
or path, pass its actual absolute reference-data URL instead.

Upload the entire generated `bikc.be` directory to the site root. This
includes `RamplotR/reference-data`, `RamplotR/app.json`, the assets in
`shinylive`, and the static HTML entry point. The data must be available
at exactly the URL passed during export. Same-origin hosting is preferred.
If serving data from another origin, configure that server's CORS policy.

`reference-index.tsv` inside each staged reference directory lists
every available group and the MD5 of its original RDS file. The browser
downloads only the files it requests and checks the downloaded bytes before
loading them. A mismatch fails closed: rebuild and redeploy the HTML app
and reference-data together. Keep the URL and deployed files in sync.

This is a **first-load** optimization, not a reduction in the number of
scientific reference datasets: users who select many datasets will download
them during their session. Ordinary Shiny and offline batch analysis still
work directly with the original checked-in local files.

## Check the export

1. Build into an empty destination or use the script to replace its app
   directory.
2. Verify `RamplotR/app.json` is substantially smaller than the original
   export and that five directories exist under `RamplotR/reference-data`.
3. Serve the output locally with `httpuv::runStaticServer("bikc.be")`;
   for a local test, export with the correct local absolute URL, not the
   production URL.
4. Check browser Network and Console tabs. Load a structure, change between
   all five reference datasets, select a per-amino-acid background, and
   exercise residue-aware and legacy classification.
5. Compare representative reference-grid MD5s and classification outputs
   between the local app and browser app. Exported source files retain
   identical bytes; independent execution in webR should still be tested.
6. Verify ordinary PNG/SVG/HTML exports and that production hosting allows
   same-origin requests to `reference-data`.

## xml2

`xml2` is not needed for routine PDB/mmCIF Ramachandran analysis.
It is needed to import official wwPDB validation XML files. A local
`shinylive::export()` warning about a missing local `xml2` is not
proof that it is absent from or supported by the selected webR package
repository. Install it locally to eliminate the local missing-package
warning. Test the XML attachment feature specifically in the exported
browser app. If no compatible webR binary exists, keep that feature in
the regular Shiny version and explain its browser limitation rather than
silently disabling it.

## First-load performance

Shinylive runs webR and Shiny entirely in the browser. A cold visitor must
download and initialize webR and the required R WebAssembly packages before
interacting with the app. The reference-data split reduces `app.json`,
but does not remove the webR boot cost.

The app no longer downloads and parses the full Plotly JavaScript bundle
before displaying the input form. Selecting **Analyse** begins fetching
Plotly while R parses the selected structure. Plot rendering waits for this
download when necessary; a prediction-only PAE heatmap can also request it.
The full bundle remains intentional: RamplotR's comparison view uses
`scattergl` for large structures, and Plotly's smaller Cartesian partial
bundle does not include this trace type.

The generated Shinylive index page uses the same small SVG favicon as the
Shiny app's embedded page.

### Measure your deployment

1. Open browser Developer Tools, Network, disable cache, and reload the
   published RamplotR page. Record the transferred size and time for
   `app.json`, `webr` assets and the downloaded `*.wasm` / R package files.
2. Without loading a structure, verify Plotly is absent in Network. Click
   **Analyse** and check that the Plotly request begins at the click.
3. Repeat with the browser cache enabled. A faster second visit points to
   downloadable webR/assets and HTTP caching as the cold-start cost. Compare
   with a locally hosted normal Shiny app if CPU startup remains slow.

On the static web server, enable Brotli or gzip for HTML, JSON, JavaScript
and other text assets, and provide sensible caching for the webR runtime,
WASM packages and unchanged reference files. Prefer `Cache-Control:
no-cache` for the generated `index.html` so deployments refresh.
Do not give mutable `app.json` or `reference-data` long immutable caching
unless the deployment uses a versioned URL, since a stale manifest paired
with new RDS files will correctly fail checksum verification.
