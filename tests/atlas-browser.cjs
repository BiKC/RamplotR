// Deterministic end-to-end smoke test of the initial Atlas inventory.
// Run after the Shiny server starts at http://127.0.0.1:8765.
// Archive data are mocked. This tests UI plumbing, not biological results.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const puppeteer = require("puppeteer-core");

(async function () {
  const chrome = [
    process.env.CHROME_BIN, "/usr/bin/google-chrome",
    "/usr/bin/google-chrome-stable", "/usr/bin/chromium"
  ].find(x => x && fs.existsSync(x));
  if (!chrome) throw new Error("Chrome/Chromium required for Atlas browser smoke test.");
  const browser = await puppeteer.launch({
    executablePath: chrome, headless: true,
    args: ["--no-sandbox", "--disable-dev-shm-usage", "--disable-gpu"]
  });
  const output = "benchmarks/output/ui-preview";
  fs.mkdirSync(output,{recursive:true});
  try {
    const page = await browser.newPage();
    await page.setViewport({width:1366,height:900});
    const requests = [];
    const errors = [];
    page.on("pageerror", error => errors.push(String(error.message || error)));
    await page.setRequestInterception(true);
    page.on("request", request => {
      const url = request.url();
      const cors = {
        "access-control-allow-origin":"*",
        "access-control-allow-methods":"GET, POST, OPTIONS",
        "access-control-allow-headers":"content-type",
        "content-type":"application/json"
      };
      if (url === "https://search.rcsb.org/rcsbsearch/v2/query") {
        if (request.method() === "OPTIONS") {
          request.respond({status:204,headers:cors,body:""});
          return;
        }
        const submitted=JSON.parse(request.postData() || "{}");
        requests.push(submitted);
        const offset=submitted.request_options.paginate.start;
        request.respond({status:200,headers:cors,body:JSON.stringify({
          total_count:73,
          result_set:[{identifier:offset===0 ? "1CRN_1" : "1UBQ_1"}]
        })});
      } else if (url === "https://www.ebi.ac.uk/pdbe/static/entry/1crn_updated.cif") {
        const mockCif = `data_1crn
#
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.entity_id
_pdbx_poly_seq_scheme.seq_id
_pdbx_poly_seq_scheme.pdb_strand_id
_pdbx_poly_seq_scheme.pdb_seq_num
_pdbx_poly_seq_scheme.pdb_ins_code
X 1 1 A 101 .
X 1 2 A 101 A
#
loop_
_pdbx_sifts_xref_db.entity_id
_pdbx_sifts_xref_db.asym_id
_pdbx_sifts_xref_db.seq_id
_pdbx_sifts_xref_db.unp_acc
_pdbx_sifts_xref_db.unp_num
_pdbx_sifts_xref_db.observed
1 X 1 P00533 50 1
1 X 2 P00533 52 1
#`;
        request.respond({status:200,
          headers:{...cors,"content-type":"text/plain"},
          body:mockCif});
      } else if (url === "https://www.ebi.ac.uk/pdbe/static/entry/1ubq_updated.cif") {
        const mockCif = `data_1ubq
#
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.entity_id
_pdbx_poly_seq_scheme.seq_id
_pdbx_poly_seq_scheme.pdb_strand_id
_pdbx_poly_seq_scheme.pdb_seq_num
_pdbx_poly_seq_scheme.pdb_ins_code
X 1 1 A 201 .
X 1 2 A 202 .
#
loop_
_pdbx_sifts_xref_db.entity_id
_pdbx_sifts_xref_db.asym_id
_pdbx_sifts_xref_db.seq_id
_pdbx_sifts_xref_db.unp_acc
_pdbx_sifts_xref_db.unp_num
_pdbx_sifts_xref_db.observed
1 X 1 P00533 51 1
1 X 2 P00533 52 1
#`;
        request.respond({status:200,
          headers:{...cors,"content-type":"text/plain"},body:mockCif});
      } else if (url === "https://data.rcsb.org/rest/v1/core/polymer_entity/1CRN/1") {
        request.respond({status:200,headers:cors,body:JSON.stringify({
          rcsb_polymer_entity_container_identifiers:{
            auth_asym_ids:["A"],
            reference_sequence_identifiers:[{
              database_name:"UniProt",database_accession:"P00533",
              reference_sequence_coverage:0.75,entity_sequence_coverage:1
            }]
          },
          rcsb_polymer_entity:{pdbx_description:"Synthetic atlas test protein"}
        })});
      } else if (url === "https://data.rcsb.org/rest/v1/core/polymer_entity/1UBQ/1") {
        request.respond({status:200,headers:cors,body:JSON.stringify({
          rcsb_polymer_entity_container_identifiers:{
            auth_asym_ids:["A"],
            reference_sequence_identifiers:[{
              database_name:"UniProt",database_accession:"P00533",
              reference_sequence_coverage:0.55,entity_sequence_coverage:0.98
            }]
          },
          rcsb_polymer_entity:{pdbx_description:"Second atlas entity"}
        })});
      } else if (url === "https://data.rcsb.org/rest/v1/core/entry/1UBQ") {
        request.respond({status:200,headers:cors,body:JSON.stringify({
          struct:{title:"Second mock entry"},
          exptl:[{method:"SOLUTION NMR"}],
          rcsb_entry_info:{resolution_combined:[]},
          rcsb_accession_info:{initial_release_date:"2001-01-01"}
        })});
      } else if (url === "https://data.rcsb.org/rest/v1/core/entry/1CRN") {
        request.respond({status:200,headers:cors,body:JSON.stringify({
          struct:{title:"Atlas mock entry"},
          exptl:[{method:"X-RAY DIFFRACTION"}],
          rcsb_entry_info:{resolution_combined:[1.5]},
          rcsb_accession_info:{initial_release_date:"1981-04-30"}
        })});
      } else {
        request.continue();
      }
    });
    await page.goto("http://127.0.0.1:8765",{
      waitUntil:"networkidle2",timeout:60000
    });
    await page.click('.nav-tabs a[data-value="atlas"]');
    await page.waitForSelector("#atlasAccession");
    await page.type("#atlasAccession","P00533");
    await page.click("#atlasDiscover");
    await page.waitForFunction(() => {
      const results=document.querySelector("#atlasResults");
      return results && results.textContent.includes("Synthetic atlas test protein") &&
        results.textContent.includes("75.0% UniProt sequence coverage");
    },{timeout:25000});
    const visible=await page.$eval("#atlasResults", el=>el.textContent);
    assert.match(visible,/73 matching entities/);
    assert.match(visible,/1 of 73 search hits retrieved/);
    assert.match(visible,/1.50 Å/);
    assert.ok(requests.length===1,"Atlas should issue one explicit search.");
    assert.equal(requests[0].query.parameters.value,"P00533");
    assert.equal(requests[0].query.parameters.operator,"exact_match");
    assert.deepEqual(requests[0].request_options.results_content_type,["experimental"]);
    assert.deepEqual(requests[0].request_options.paginate,{start:0,rows:50});
    assert.deepEqual(requests[0].request_options.sort,
      [{sort_by:"rcsb_id",direction:"asc"}]);
    // A selected experimental entity can verify exact SIFTS positions
    // without guessing interior numbering or losing insertion codes.
    await page.click(".ram-atlas-verify");
    await page.waitForFunction(() => {
      const results=document.querySelector("#atlasResults");
      return results && results.textContent.includes("Verified exact SIFTS:") &&
        results.textContent.includes("2 distinct PDB residues") &&
        results.textContent.includes("2 observed");
    },{timeout:25000});
    await page.click("#atlasLoadMore");
    await page.waitForFunction(() => {
      const results=document.querySelector("#atlasResults");
      return results && results.textContent.includes("Second atlas entity") &&
        results.textContent.includes("2 of 73 search hits retrieved");
    },{timeout:25000});
    assert.equal(requests.length,2,"Load more should issue exactly one new page request.");
    assert.deepEqual(requests[1].request_options.paginate,{start:1,rows:50});
    const afterPaging=await page.$eval("#atlasResults",el=>el.textContent);
    assert.match(afterPaging,/2 enriched polymer entities in 2 PDB entries/);
    assert.match(afterPaging,/55.0% UniProt sequence coverage/);
    const verifier=await page.$$(".ram-atlas-verify");
    assert.equal(verifier.length,2,"Expected one verifier per experimental entity.");
    await verifier[1].click();
    await page.waitForFunction(() => {
      const text=document.querySelector("#atlasResults")?.textContent || "";
      return text.includes("2 verified experimental entities") &&
        text.includes("3 distinct observed UniProt positions") &&
        text.includes("Export verified residue mapping CSV");
    },{timeout:25000});
    assert.ok(await page.$("#downloadAtlasCanonical"));
    assert.ok(await page.$("#downloadAtlasSupport"));
    await page.screenshot({path:path.join(output,"atlas-inventory-desktop.png"),
      fullPage:true});
    await page.setViewport({width:390,height:844});
    await page.screenshot({path:path.join(output,"atlas-inventory-mobile.png"),
      fullPage:true});
    // Results are not permitted to act as conformation-cluster labels.
    assert.doesNotMatch(visible,/experimental state [1-9]/i);
    assert.ok(!errors.length,"Atlas browser JavaScript errors: "+errors.join(" | "));
    console.log("Atlas inventory browser smoke test passed.");
  } finally {
    await browser.close();
  }
})().catch(error => {console.error(error);process.exitCode=1;});
