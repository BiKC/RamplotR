# Atlas browser archive connectivity

RamplotR is deployed as static Shinylive at
[https://bikc.be/RamplotR/](https://bikc.be/RamplotR/).
Archive queries originate in the **visitor's browser**: a command-line
`curl`, a local Shiny test or GitHub-hosted Node HTTP request does not
prove that the site is permitted to fetch RCSB/PDBe under its browser
security policy.

## In-app check

Go to **Atlas → Check browser access to RCSB and PDBe → Test archive
connections**. The check is opt-in, leaves the current analysis untouched
and tests three concrete endpoints with the same parsing helpers as Atlas:

| Step | Archive URL | Success criterion |
| --- | --- | --- |
| Experimental search | `POST https://search.rcsb.org/rcsbsearch/v2/query` | Valid JSON search response for UniProt P69441 (using the Atlas experimental-only query builder). |
| Polymer-entity metadata | `GET https://data.rcsb.org/rest/v1/core/polymer_entity/4AKE/1` | Valid RCSB polymer metadata JSON. |
| Exact experimental SIFTS | `GET https://www.ebi.ac.uk/pdbe/static/entry/4ake_updated.cif` | Parse exact SIFTS P69441 residue mappings, observed records and first-model Cα coordinates, using the real Atlas browser parser. |

The read-only probe uses `credentials: "omit"`, disables response caching,
imposes request timeouts and displays every failure independently.
Successful checks do **not** guarantee every archive entry or larger
coordinate file is accessible.

The browser's network error alone cannot reliably distinguish blocked CORS
from network interruption, content-security policy, an extension, TLS failure
or offline mode. Accordingly, this is reported as **network or CORS** rather
than claiming a specific cause.

## Troubleshooting

- **RCSB Search fails, RCSB metadata works:** check the browser developer
  Network console for the `OPTIONS` preflight and `POST`, including
  `Access-Control-Allow-Origin` and allowed methods/headers.
- **RCSB metadata fails:** inspect the GET response status and content;
  check whether the browser, institutional network, or privacy extension
  blocks `data.rcsb.org`.
- **PDBe download fails:** check `www.ebi.ac.uk` response status, CORS
  headers and file size. A successful HTTP 200 with no parseable exact
  SIFTS or model-1 Cα records is **not** a passing check.
- **Works locally but not at bikc.be:** inspect the page's actual origin,
  browser console, iframe and content-security policy; repeat from the
  deployed site. Static hosting can have different restrictions than a
  localhost Shiny session.
- **Timeout/temporary 5xx:** retry when the provider is available.
  Never silently replace verified residue mapping with guessed numbering
  or route personally uploaded protein data to an unrelated proxy.

## Automated verification

`.github/workflows/atlas-browser-origin.yml` is triggered **once when a
pull request becomes Ready for review**, with a manual fallback. It visits
the real `https://bikc.be/RamplotR/` in headless Chromium, injects the
checked-out branch's actual Atlas parsers and executes the three probes from
that origin. This tests browser-origin archive access **even before the
latest application UI has been deployed**. Results, elapsed times and
per-endpoint diagnostics are saved in
`ramplotr-public-browser-archive-check` artifacts.

This is distinct from the local end-to-end browser test, which uses mocked
archive responses to deterministically test UI state and failure handling.
The automated public-origin check cannot prove a full newly deployed
Shinylive workflow works: that still requires publishing the latest build
and using its UI in an actual browser.

Scientific methods and thresholds are unaffected by these network checks.
