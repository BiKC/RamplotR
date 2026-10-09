// Deterministic end-to-end smoke test of the initial Atlas inventory.
// Run after the Shiny server starts at http://127.0.0.1:8765.
// Archive data are mocked. This tests UI plumbing, not biological results.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const puppeteer = require("puppeteer-core");

// Realistic synthetic common cores exercise browser mmCIF C-alpha extraction
// and the Shiny distance-map grouping without depending on live PDBe APIs.
function syntheticUpdatedCif(pdbOffset,shifted) {
  const residues=Array.from({length:40},(_,i)=>i+1);
  const sequence=residues.map(i=>`X 1 ${i} A ${pdbOffset+i-1} . ALA`).join("\n");
  const sifts=residues.map(i=>`1 X ${i} P00533 ${i+49} 1`).join("\n");
  const atoms=residues.flatMap(i=>{
    const base=[1.5*i,0.5*Math.sin(i),0.5*Math.cos(i)];
    const p=[base[0]+(shifted && i>21 ? 8 : 0),base[1],base[2]];
    const N=[p[0]-0.2,p[1]+0.3,p[2]+0.1];
    const C=[p[0]+0.2,p[1]-0.1,p[2]+0.25];
    if(shifted && i>=15 && i<=19) {
      N[1]+=0.45; N[2]-=0.1;
      C[1]-=0.3; C[2]+=0.3;
    }
    return [["N",N],["CA",p],["C",C]].map(([name,xyz])=>
      "ATOM "+name+" X "+i+" "+xyz.map(x=>x.toFixed(4)).join(" ")+" 1 .");
  }).join("\n");
  return `data_test
#
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.entity_id
_pdbx_poly_seq_scheme.seq_id
_pdbx_poly_seq_scheme.pdb_strand_id
_pdbx_poly_seq_scheme.pdb_seq_num
_pdbx_poly_seq_scheme.pdb_ins_code
_pdbx_poly_seq_scheme.mon_id
${sequence}
#
loop_
_pdbx_sifts_xref_db.entity_id
_pdbx_sifts_xref_db.asym_id
_pdbx_sifts_xref_db.seq_id
_pdbx_sifts_xref_db.unp_acc
_pdbx_sifts_xref_db.unp_num
_pdbx_sifts_xref_db.observed
${sifts}
#
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
_atom_site.label_alt_id
${atoms}
#`;
}
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
        if(submitted.query?.parameters?.value==="P69441") {
          request.respond({status:200,headers:cors,
            body:JSON.stringify({total_count:1,
              result_set:[{identifier:"4AKE_1"}]})});
          return;
        }
        requests.push(submitted);
        const offset=submitted.request_options.paginate.start;
        request.respond({status:200,headers:cors,body:JSON.stringify({
          total_count:73,
          result_set:[{identifier:offset===0 ? "1CRN_1" : "1UBQ_1"}]
        })});
      } else if (url === "https://data.rcsb.org/rest/v1/core/polymer_entity/4AKE/1") {
        request.respond({status:200,headers:cors,body:JSON.stringify({
          rcsb_polymer_entity_container_identifiers:{
            entry_id:"4AKE",entity_id:"1",auth_asym_ids:["A"]
          }
        })});
      } else if (url === "https://www.ebi.ac.uk/pdbe/static/entry/4ake_updated.cif") {
        request.respond({status:200,
          headers:{...cors,"content-type":"text/plain"},
          body:syntheticUpdatedCif(101,false).replace(/P00533/g,"P69441")});
      } else if (url === "https://www.ebi.ac.uk/pdbe/static/entry/1crn_updated.cif") {
        const mockCif = syntheticUpdatedCif(101,false);
        request.respond({status:200,
          headers:{...cors,"content-type":"text/plain"},
          body:mockCif});
      } else if (url === "https://www.ebi.ac.uk/pdbe/static/entry/1ubq_updated.cif") {
        const mockCif = syntheticUpdatedCif(201,true);
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
    // The network diagnostic runs only when explicitly opened/clicked and
    // must not mutate or trigger the Atlas search itself.
    await page.click(".ram-atlas-network-check > summary");
    await page.click("#atlasConnectivityRun");
    await page.waitForFunction(()=>{
      const panel=document.querySelector("#atlasConnectivityStatus");
      return panel && panel.textContent.includes("All 3 endpoints passed");
    },{timeout:25000});
    const probe=await page.$eval("#atlasConnectivityStatus",
      el=>el.textContent);
    assert.match(probe,/Browser origin: http:\/\/127\.0\.0\.1:8765/);
    assert.match(probe,/PDBe updated mmCIF/);
    assert.equal(await page.$eval(
      ".ram-atlas-connectivity-list li",nodes=>nodes.length),3);
    assert.equal(requests.length,0,
      "Connectivity test should not trigger a real Atlas cohort search.");
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
        results.textContent.includes("40 distinct PDB residues") &&
        results.textContent.includes("40 observed");
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
        text.includes("40 distinct observed UniProt positions") &&
        text.includes("Export verified residue mapping CSV");
    },{timeout:25000});
    assert.ok(await page.$("#downloadAtlasCanonical"));
    assert.ok(await page.$("#downloadAtlasSupport"));
    await page.waitForSelector("#atlasGeometryEntities",{timeout:15000});
    await page.waitForFunction(() => {
      const review=document.querySelector("#atlasConstructReview");
      const data=document.querySelector("#atlasConstructTable");
      return review && review.textContent.includes("Construct and sequence-chemistry review") &&
        data && data.textContent.includes("Chemistry checked");
    },{timeout:20000});
    const reviewText=await page.$eval("#atlasConstructReview",el=>el.textContent);
    assert.match(reviewText,/2 experimental-entity pairs|1 experimental-entity pair/i);
    assert.ok(await page.$("#downloadAtlasConstruct"));
    await page.click("#atlasRunGeometry");
    await page.waitForFunction(() => {
      const panel=document.querySelector("#atlasGeometrySummary");
      return panel && panel.textContent.includes(
        "2 exploratory geometric group(s)") &&
        panel.textContent.includes("40 common observed UniProt");
    },{timeout:25000});
    const geometryTable=await page.$eval("#atlasGeometryTable",
      el=>el.textContent);
    assert.match(geometryTable,/1CRN_1/);
    assert.match(geometryTable,/1UBQ_1/);
    try {
      await page.waitForFunction(() => {
        const plot=document.querySelector("#atlasGeometryPlot");
        const img=plot && plot.querySelector("img");
        return !!(plot && (plot.querySelector("canvas") ||
          (img && img.complete && img.naturalWidth>0)));
      },{timeout:25000});
    } catch(e) {
      console.error("Atlas geometry plot diagnostics:",
        JSON.stringify(await page.evaluate(() => ({
          plot:document.querySelector("#atlasGeometryPlot")?.outerHTML?.slice(0,1500),
          result:document.querySelector("#atlasGeometrySummary")?.textContent,
          notices:[...document.querySelectorAll(".shiny-notification")]
            .map(node=>node.textContent)
        }))));
      throw e;
    }

    try {
      await page.waitForSelector("#atlasSwitchThreshold",{timeout:15000});
    } catch(error) {
      console.error("Atlas switch panel diagnostics:",JSON.stringify(
        await page.evaluate(()=>({
          geometry:document.querySelector("#atlasGeometrySummary")?.textContent,
          switchPanel:document.querySelector("#atlasSwitchPanel")?.outerHTML?.slice(0,2000),
          errors:[...document.querySelectorAll(".shiny-output-error")]
            .map(node=>node.textContent?.slice(0,400)),
          messages:[...document.querySelectorAll(".shiny-notification")]
            .map(node=>node.textContent)
        }))));
      throw error;
    }
    await page.click("#atlasRunSwitch");
    await page.waitForFunction(() => {
      const panel=document.querySelector("#atlasSwitchSummary");
      return panel && panel.textContent.includes(
        "canonical positions have both valid") &&
        panel.textContent.includes("candidate residue");
    },{timeout:25000});
    const switchText=await page.$eval("#atlasSwitchSummary",el=>el.textContent);
    assert.match(switchText,/1CRN_1 versus 1UBQ_1/);
    assert.ok(await page.$("#downloadAtlasSwitch"));
    await page.waitForFunction(() => {
      const plot=document.querySelector("#atlasSwitchPlot img");
      return plot && plot.complete && plot.naturalWidth>0;
    },{timeout:25000});
    // The Atlas table must lead to a residue inspector before any model
    // navigation occurs. The synthetic mappings deliberately differ from
    // real PDB author numbering, so the test does not claim paired 3D focus.
    await page.waitForSelector("#atlasSwitchResidues tbody tr",{timeout:25000});
    await page.click("#atlasSwitchResidues tbody tr:first-child");
    await page.waitForFunction(() => {
      const inspector=document.querySelector("#atlasSwitchResidueInspector");
      return inspector && inspector.textContent.includes("UniProt residue") &&
        inspector.textContent.includes("Chain A, residue");
    },{timeout:15000});
    const inspectorText=await page.$eval("#atlasSwitchResidueInspector",
      el=>el.textContent);
    assert.match(inspectorText,/Inspect both structures in 2D\/3D/);
    assert.match(inspectorText,/Chain A, residue/);
    assert.match(inspectorText,/first-model experimental torsions/);
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
