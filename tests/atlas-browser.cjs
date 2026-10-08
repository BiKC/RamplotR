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
        requests.push(JSON.parse(request.postData() || "{}"));
        request.respond({status:200,headers:cors,body:JSON.stringify({
          total_count:73,
          result_set:[{identifier:"1CRN_1"}]
        })});
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
    assert.match(visible,/first 1 were requested/);
    assert.match(visible,/1.50 Å/);
    assert.ok(requests.length===1,"Atlas should issue one explicit search.");
    assert.equal(requests[0].query.parameters.value,"P00533");
    assert.equal(requests[0].query.parameters.operator,"exact_match");
    assert.deepEqual(requests[0].request_options.results_content_type,["experimental"]);
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
