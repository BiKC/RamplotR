# Thin Shinylive deployment

The normal app and batch CLI continue using all five local RDS reference
directories. A thin browser deployment keeps the same reference options,
but moves the original, unchanged RDS files out of the initial app.json.

## Build

From the repository root, install `shinylive` and run the optimized export:

```bash
Rscript scripts/export-shinylive.R bikc.be https://bikc.be/RamplotR/reference-data
```

On Windows, use `Rscript.exe` from your R installation if it is not on
`PATH`. The second argument must be the public URL of the deployed
reference-data directory. The ordinary `shinylive::export()` command
still creates a full package containing every bundled RDS file, so use the
repository script for the smaller browser deployment.

Upload the **contents** of the generated `bikc.be` directory to the site document root. This
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

## one.com hosting

one.com supports Apache `.htaccess` files but restricts some directives.
The export includes three optional, scoped configurations. They **do not**
change the website root `.htaccess` or the setup of other apps:

| Generated file | Effect |
| --- | --- |
| `RamplotR/.htaccess` | Enables Brotli or gzip for the browser app's HTML, JSON, JS and CSS if the matching Apache module is available. Revalidates `index.html` and `app.json`; caches local assets for one day. |
| `RamplotR/reference-data/.htaccess` | Revalidates the RDS files between deployments. They are already gzip-compressed R objects and should not be recompressed. |
| `shinylive/.htaccess` | Optionally compresses shared webR and WASM assets and caches static assets for one day, but revalidates metadata. An existing shared `.htaccess` is **not overwritten**. If one exists, review and merge the template under `config/onecom/` manually. |

The files are in `config/onecom/` if you need to inspect or adjust them.
One.com may not allow every `mod_brotli`, `mod_deflate` or `mod_headers`
directive. Conditional module blocks mean missing modules are skipped, but
they do not bypass one.com's hosting restrictions. If you get a 500 error
after uploading, remove the generated `.htaccess` files and ask one.com
support whether those directives are permitted on your hosting plan.
Your Shiny app does not require these rules to function.

Do not enable the WordPress-specific Performance Cache plugin as a
requirement for this static app. It is separate from ordinary HTTP caching.

### Confirm HTTP compression and caching

In a terminal, inspect response headers (use `curl.exe` on Windows if
PowerShell's `curl` alias is active):

```bash
curl -I -H "Accept-Encoding: br,gzip" https://bikc.be/RamplotR/app.json
curl -I -H "Accept-Encoding: br,gzip" https://bikc.be/RamplotR/favicon.svg
```

The `app.json` request should show `Cache-Control: no-cache` when `mod_headers`
is enabled. For text resources, `Content-Encoding: br` or `gzip` indicates
compression is active. No such header means the host is not compressing
that response or the requested file is not found; check with the browser
Network panel using the correct asset URL from the generated page.
For WASM and other webR assets, inspect the actual paths reported in the
Network panel rather than assuming their locations.

To assess a first visit, disable browser cache and record network transfers
separately from webR startup and R-package initialization. Then reload
with caching enabled. A long download points to hosting and network costs;
slow initialization after all downloads points to webR/package startup
or device CPU. This is useful before making further application changes.
