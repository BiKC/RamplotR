// End-to-end predicted-model smoke tests against the running Shiny app.
// The B factors below are synthetic test data, not a real ESMFold prediction.
const fs = require("node:fs");
const path = require("node:path");
const assert = require("node:assert/strict");
const puppeteer = require("puppeteer-core");
(async function () {
  const output = "benchmarks/output/ui-preview";
  fs.mkdirSync(output, {recursive:true});
  const raw = fs.readFileSync(path.join(output,"1CRN.pdb"),"utf8");
  const pdb = raw.split(/\r?\n/).map(line => line.startsWith("ATOM  ")
    ? line.slice(0,60) + " 87.00" + line.slice(66) : line).join("\n");
  const fixture = path.resolve(output,"synthetic-esmfold.pdb");
  fs.writeFileSync(fixture,pdb);
  const pdbSeed2 = raw.split(/\r?\n/).map(line => line.startsWith("ATOM  ")
    ? line.slice(0,60) + " 67.00" + line.slice(66) : line).join("\n");
  const fixtureSeed2 = path.resolve(output,"synthetic-esmfold-seed2.pdb");
  fs.writeFileSync(fixtureSeed2,pdbSeed2);
  const atomLines = pdb.split("\n").filter(line => line.startsWith("ATOM  "));
  const unique = [...new Set(atomLines.map(line =>
    line.slice(21,22) + "|" + line.slice(22,26).trim() + "|" +
    line.slice(26,27).trim()))];
  const n = unique.length;
  const pae = Array.from({length:n}, (_,i) =>
    Array.from({length:n},(_,j) => Math.abs(i-j)*0.5));
  const af2 = path.resolve(output,"synthetic-af2-pae.json");
  const af3 = path.resolve(output,"synthetic-af3-confidence.json");
  const af3Summary = path.resolve(output,"synthetic-af3-summary.json");
  fs.writeFileSync(af2,JSON.stringify([{
    predicted_aligned_error:pae,max_predicted_aligned_error:31.75
  }]));
  fs.writeFileSync(af3,JSON.stringify({
    pae, token_chain_ids:unique.map(key=>key.split("|")[0]),
    token_res_ids:unique.map(key=>Number(key.split("|")[1])),
    atom_chain_ids:atomLines.map(line=>line.slice(21,22)),
    atom_plddts:atomLines.map(()=>88)
  }));
  fs.writeFileSync(af3Summary,JSON.stringify({ptm:0.82,iptm:0.71}));
  const chrome = [process.env.CHROME_BIN,"/usr/bin/google-chrome",
    "/usr/bin/google-chrome-stable","/usr/bin/chromium"].find(x=>x&&fs.existsSync(x));
  if (!chrome) throw Error("Chrome or Chromium is required.");
  const browser = await puppeteer.launch({
    executablePath:chrome,headless:true,
    args:["--no-sandbox","--disable-dev-shm-usage","--disable-gpu"]
  });
  try {
    const page = await browser.newPage();
    page.on("pageerror",error=>console.error("Browser error:",error.message));
    await page.setViewport({width:1366,height:900,deviceScaleFactor:1});
    await page.goto("http://127.0.0.1:8765",{
      waitUntil:"networkidle2",timeout:60000});
    await page.evaluate(() => document.querySelector(
      'input[name="inputSource"][value="upload"]').click());
    await page.waitForFunction(()=>!document.getElementById("ram-upload-wrap")
      .classList.contains("is-hidden"));
    const upload=await page.$("#structfile");
    await upload.uploadFile(fixture);
    // Select prediction provenance through the same visible native dropdown a
    // researcher uses. A Selectize-backed hidden input may silently revert
    // scripted changes and analyse an AF/ESM prediction as experimental.
    async function chooseSource(value) {
      const settings=await page.$("#ram-prediction-upload");
      const opened=await page.evaluate(el=>el.open,settings);
      if (!opened) await page.click("#ram-prediction-upload > summary");
      await page.select("#predictionSource",value);
      await page.waitForFunction(expected =>
        window.Shiny && window.Shiny.shinyapp &&
        window.Shiny.shinyapp.$inputValues.predictionSource === expected,
        {timeout:10000},value);
    }
    await chooseSource("esmfold");
    await page.waitForFunction(()=>document.getElementById("ram-confidence-sidecars")
      .classList.contains("is-hidden"));
    await new Promise(done=>setTimeout(done,1300));
    await page.click("#submit");
    try {
      await page.waitForFunction(() => {
        const panel=document.querySelector("#ram-confidence-panel");
        return panel && panel.textContent.includes("ESMFold") &&
          panel.textContent.includes("87.0") &&
          document.querySelectorAll(".ram-confidence-mini-cell").length>0;
      },{timeout:25000});
    } catch (error) {
      const diagnostics=await page.evaluate(() => ({
        source: document.querySelector('input[name="inputSource"]:checked')?.value,
        predictionSource: document.querySelector("#predictionSource")?.value,
        predictionSettingsOpen: document.querySelector("#ram-prediction-upload")?.open,
        fileName:document.querySelector("#structfile")?.files?.[0]?.name,
        fileInput: window.Shiny?.shinyapp?.$inputValues?.["structfile:shiny.file"],
        shinyInput: window.Shiny?.shinyapp?.$inputValues?.predictionSource,
        panel:document.querySelector("#predictionPanel")?.textContent?.slice(0,1200),
        plot:document.querySelector("#ram-current-structure")?.textContent,
        notifications:[...document.querySelectorAll(".shiny-notification")]
          .map(el=>el.textContent)
      }));
      console.error("Prediction upload diagnostics:",JSON.stringify(diagnostics));
      await page.screenshot({path:path.join(output,"prediction-debug.png"),fullPage:true});
      throw error;
    }
    await page.waitForSelector("#ram-confidence-panel > summary",{timeout:12000});
    const confidencePanelOpen=await page.$eval("#ram-confidence-panel",
      panel=>panel.open);
    if(!confidencePanelOpen)
      await page.click("#ram-confidence-panel > summary");
    await page.waitForSelector("#counterpartAccession",
      {visible:true,timeout:12000});
    const counterpartUi=await page.evaluate(() => ({
      accession:document.getElementById("counterpartAccession")?.value || "",
      lookup:Boolean(document.getElementById("lookupCounterparts")),
      use:Boolean(document.getElementById("useCounterpart"))
    }));
    assert.equal(counterpartUi.accession,"",
      "Uploaded predictions should not invent a UniProt accession.");
    assert.ok(counterpartUi.lookup && counterpartUi.use,
      "Prediction confidence should expose experimental counterpart controls.");
    await page.$eval("#counterpartAccession",input=>{
      input.value="not an accession";
      input.dispatchEvent(new Event("input",{bubbles:true}));
      input.dispatchEvent(new Event("change",{bubbles:true}));
    });
    await page.click("#lookupCounterparts");
    await page.waitForFunction(() =>
      [...document.querySelectorAll(".shiny-notification")]
        .some(node=>node.textContent.includes("valid UniProt accession")),
      {timeout:12000});

    assert.equal(await page.$("#ram-pae-plot"),null,
      "ESMFold without PAE must not invent an error map");
    await page.screenshot({path:path.join(output,"prediction-esmfold.png"),
      fullPage:true});
    await page.click("#ram-sequence-panel > summary");
    await page.waitForFunction(() => {
      const buttons = Array.from(document.querySelectorAll("#sequenceView .ram-seq-res"));
      return buttons.length > 20 &&
        buttons.some(button => button.dataset.plddt === "87.0") &&
        buttons.some(button => button.querySelector(".ram-seq-plddt")?.textContent === "87") &&
        document.querySelector(".ram-sequence-confidence-key");
    },{timeout:20000});
    const confidence = await page.$eval("#sequenceView .ram-seq-res[data-plddt]",
      button => ({
        label:button.querySelector(".ram-seq-plddt").textContent,
        underline:getComputedStyle(button).boxShadow,
        tooltip:button.title
      }));
    assert.equal(confidence.label,"87");
    assert.ok(confidence.tooltip.includes("pLDDT 87.0"));
    assert.ok(confidence.underline !== "none",
      "Each confidence-bearing residue must have a separate colour strip");
    await page.screenshot({
      path:path.join(output,"prediction-sequence-confidence.png"),fullPage:true
    });
    await page.click("#ram-sequence-panel > summary");

    // Independently verify AF2 monomer PAE and UI linkage.
    await chooseSource("alphafold2");
    await page.waitForFunction(()=>!document.getElementById("ram-confidence-sidecars")
      .classList.contains("is-hidden"));
    await (await page.$("#predictionJson")).uploadFile(af2);
    await new Promise(done=>setTimeout(done,1200));
    await page.click("#submit");
    await page.waitForFunction(()=> {
      const el=document.querySelector("#ram-confidence-panel");
      return el && el.textContent.includes("AlphaFold");
    },{timeout:90000});
    await page.click("#ram-confidence-panel summary");
    await page.waitForFunction(()=>{
      const plot=document.querySelector("#ram-pae-plot");
      return plot && plot.data && plot.data.length>=2 &&
        plot.data[0].z.length>2;
    },{timeout:30000});
    const map=await page.evaluate(() => {
      const plot=document.querySelector("#ram-pae-plot");
      return {n:plot.data[0].z.length,first:plot.data[0].z[0][1]};
    });
    assert.equal(map.n,n);
    assert.equal(map.first,0.5);
    await page.screenshot({path:path.join(output,"prediction-af2-pae.png"),
      fullPage:true});

    // Explicit AF3 upload uses per-atom confidence and per-token PAE. Both
    // must map to the same matching chain/residue, without reusing an AF2 label.
    await chooseSource("alphafold3");
    await (await page.$("#predictionJson")).uploadFile(af3);
    await (await page.$("#predictionSummaryJson")).uploadFile(af3Summary);
    await new Promise(done=>setTimeout(done,1300));
    await page.click("#submit");
    await page.waitForFunction(()=>{
      const el=document.querySelector("#ram-confidence-panel");
      return el && el.textContent.includes("AlphaFold 3") &&
        el.textContent.includes("88.0");
    },{timeout:90000});
    await page.click("#ram-confidence-panel summary");
    await page.waitForFunction(()=>{
      const p=document.getElementById("ram-pae-plot");
      return p && p.data && p.data[0].z.length>2;
    },{timeout:30000});
    await page.waitForFunction(()=>document.querySelector(
      ".ram-confidence-metrics")?.textContent.includes("0.82"));
    await page.screenshot({path:path.join(output,"prediction-af3.png"),
      fullPage:true});

    // Prediction ensemble: two separate model files, deliberately carrying
    // different synthetic pLDDT values. This tests multi-file upload,
    // aggregation, variability rendering and map -> residue selection.
    await page.click('.nav-tabs a[data-value="summary"]');
    await page.waitForSelector("#ram-ensemble-panel",{timeout:20000});
    const ensemblePanel = await page.$("#ram-ensemble-panel");
    if (!(await page.evaluate(el=>el.open,ensemblePanel)))
      await page.click("#ram-ensemble-panel > summary");
    await page.waitForSelector("#predictionEnsembleSource",{timeout:12000});
    await page.select("#predictionEnsembleSource","esmfold");
    const ensembleUpload = await page.$("#predictionEnsembleFiles");
    await ensembleUpload.uploadFile(fixture,fixtureSeed2);
    await page.waitForFunction(() => {
      const input=document.getElementById("predictionEnsembleFiles");
      return input && input.files && input.files.length===2;
    },{timeout:10000});
    await new Promise(done=>setTimeout(done,1200));
    await page.click("#calculatePredictionEnsemble");
    await page.waitForFunction(() => {
      const panel=document.querySelector("#predictionEnsembleSummary");
      return panel && panel.textContent.includes("2 models analysed") &&
        document.querySelectorAll(".ram-ensemble-cell").length>20 &&
        document.querySelectorAll("#predictionEnsembleRows tbody tr").length>0;
    },{timeout:35000});
    const ensembleState=await page.evaluate(() => {
      const summary=document.querySelector("#predictionEnsembleSummary").textContent;
      const rows=Array.from(document.querySelectorAll(
        "#predictionEnsembleRows tbody tr"));
      const first=rows[0] ? Array.from(rows[0].querySelectorAll("td"))
        .map(td=>td.textContent.trim()) : [];
      return {
        summary,
        cells:document.querySelectorAll(".ram-ensemble-cell").length,
        first
      };
    });
    const disagreement = ensembleState.summary.match(
      /(\d+) residues with pLDDT SD ≥10/);
    assert.ok(disagreement && Number(disagreement[1]) > 0,
      "Synthetic 20-point pLDDT differences must produce a nonzero disagreement count.");
    assert.ok(ensembleState.cells>20,
      "Prediction ensemble variability map should cover the protein.");
    await page.click(".ram-ensemble-cell");
    await page.waitForFunction(() =>
      document.querySelector("#selectedResidueInfo strong"),
      {timeout:15000});
    await page.screenshot({
      path:path.join(output,"prediction-ensemble.png"),fullPage:true
    });

    console.log("Live ESMFold, AlphaFold2, AlphaFold3 and prediction-ensemble views passed.");
  } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
