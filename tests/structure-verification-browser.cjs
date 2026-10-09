// Experimental verification live Shiny smoke test. Fixtures are synthetic, not claims about
// experimental map fit or official wwPDB assessments.
const fs=require("node:fs"),path=require("node:path");
const assert=require("node:assert/strict");
const puppeteer=require("puppeteer-core");
(async function(){
  const output=path.resolve("benchmarks/output/ui-preview");
  const pdb=path.join(output,"1CRN.pdb"),nmr=path.join(output,"1D3Z.pdb");
  const atom=fs.readFileSync(pdb,"utf8").split(/\r?\n/)
    .find(line=>line.startsWith("ATOM  "));
  const chain=atom.slice(21,22),resno=atom.slice(22,26).trim();
  const resn=atom.slice(17,20).trim(),insert=atom.slice(26,27).trim();
  const xml=path.join(output,"synthetic-official-geometry.xml");
  fs.writeFileSync(xml,[
    '<?xml version="1.0"?>',
    '<wwPDB-validation-information>',
    '<ModelledSubgroup model="1" chain="'+chain+'" resnum="'+resno+
      '" icode="'+insert+'" resname="'+resn+
      '" rama="Allowed" rota="OUTLIER">',
    '<clash atom="CB" overlap="0.5"/>',
    '</ModelledSubgroup></wwPDB-validation-information>'
  ].join("\n"));
  const map=path.join(output,"synthetic-display-only.map");
  fs.writeFileSync(map,Buffer.alloc(2048));
  const chrome=[process.env.CHROME_BIN,"/usr/bin/google-chrome",
    "/usr/bin/google-chrome-stable","/usr/bin/chromium"]
    .find(f=>f&&fs.existsSync(f));
  if(!chrome)throw Error("Chrome/Chromium not found");
  const browser=await puppeteer.launch({executablePath:chrome,headless:true,
    args:["--no-sandbox","--disable-dev-shm-usage","--disable-gpu"]});
  try{
    const page=await browser.newPage();
    page.on("pageerror",error=>console.error("Page error:",error.stack||error));
    await page.setViewport({width:1366,height:900});
    await page.goto("http://127.0.0.1:8765",
                    {waitUntil:"networkidle2",timeout:60000});
    await page.evaluate(()=>document.querySelector(
      'input[name="inputSource"][value="upload"]').click());
    await page.waitForFunction(()=>
      !document.getElementById("ram-upload-wrap").classList.contains("is-hidden"));
    await (await page.$("#structfile")).uploadFile(pdb);
    await new Promise(resolve=>setTimeout(resolve,1400));
    await page.click("#submit");
    await page.waitForFunction(()=>{
      const node=document.querySelector("#ram-geometry-panel");
      return node&&node.querySelector("summary")&&
        node.textContent.includes("Extended structure verification");
    },{timeout:90000});
    // Opening a native details element directly avoids intermittent
    // Chromium clickable-point failures after async Shiny layout changes.
    await page.$eval("#ram-geometry-panel",node=>{ node.open=true; });
    await page.waitForSelector("#validationXml");
    await (await page.$("#validationXml")).uploadFile(xml);
    await new Promise(resolve=>setTimeout(resolve,1400));
    await page.click("#confirmValidationSource");
    await page.waitForFunction(() =>
      window.Shiny && window.Shiny.shinyapp &&
      window.Shiny.shinyapp.$inputValues &&
      window.Shiny.shinyapp.$inputValues.confirmValidationSource === true,
      {timeout:15000});
    await page.click("#attachValidation");
    await page.waitForFunction(()=>
      document.querySelector(".ram-official-summary")&&
      document.querySelector(".ram-official-summary").textContent
        .includes("Matched 1 of"),{timeout:30000});
    await page.screenshot({path:path.join(output,"verification-wwpdb.png"),
                           fullPage:true});
    // Mock only volume parsing to test local map controls deterministically.
    // Physical map alignment still requires a genuine CCP4 map and review.
    await page.evaluate(()=>{
      const stage=window.getNGLStage&&window.getNGLStage("NGL");
      if(!stage)throw Error("NGL stage unavailable.");
      window.__ramDensityChecks={loaded:0,levels:[],removed:0};
      const originalLoad=stage.loadFile.bind(stage);
      const originalRemove=stage.removeComponent.bind(stage);
      stage.__ramOriginalLoadFile=originalLoad;
      stage.__ramOriginalRemoveComponent=originalRemove;
      stage.loadFile=async(file,opts)=>{
        // Only simulate volume parsing. NGLVieweR also uses loadFile for
        // ordinary structure updates and must not be intercepted.
        if(opts&&opts.ext==="ccp4") {
          window.__ramDensityChecks.loaded++;
          return {addRepresentation:(name,params)=>({
            setParameters:values=>window.__ramDensityChecks.levels.push(values)
          })};
        }
        return originalLoad(file,opts);
      };
      stage.removeComponent=comp=>{
        if(comp && typeof comp.addRepresentation==="function" &&
           !comp.structure && !comp.volume) {
          window.__ramDensityChecks.removed++;
          return;
        }
        return originalRemove(comp);
      };
    });
    await page.click(".ram-density-panel summary");
    await page.waitForSelector("#ram-density-file",{timeout:15000});
    const browserLocal=await page.$eval("#ram-density-file",input=>
      !input.classList.contains("shiny-bound-input") &&
      !input.closest(".shiny-input-container"));
    assert.equal(browserLocal,true,
      "CCP4/MRC map files must remain in the browser and never enter Shiny's upload bindings.");
    await (await page.$("#ram-density-file")).uploadFile(map);
    await page.waitForFunction(()=>{
      const picker=document.getElementById("ram-density-file");
      return picker&&picker.files&&picker.files.length===1 &&
        document.getElementById("ram-density-status").textContent.includes("Selected ");
    });
    await page.click("#ram-density-load");
    try {
      await page.waitForFunction(()=>document.getElementById("ram-density-status")
        .textContent.includes("displayed at"),{timeout:12000});
    } catch(error) {
      const diagnostics=await page.evaluate(()=>({
        status:document.getElementById("ram-density-status").textContent,
        filename:document.getElementById("ram-density-file").files[0]&&
                 document.getElementById("ram-density-file").files[0].name,
        counters:window.__ramDensityChecks,
        stageAvailable:typeof window.getNGLStage==="function" &&
                       !!window.getNGLStage("NGL")
      }));
      console.error("Experimental verification map controls:",JSON.stringify(diagnostics));
      await page.screenshot({path:path.join(output,"verification-map-failure.png"),
                             fullPage:true});
      throw error;
    }
    await page.$eval("#ram-density-level",field=>{
      field.value="3.25";field.dispatchEvent(new Event("change",{bubbles:true}));
    });
    const overlay=await page.evaluate(()=>window.__ramDensityChecks);
    assert.equal(overlay.loaded,1);
    assert.equal(overlay.levels[0].isolevel,3.25);
    await page.click("#ram-density-clear");
    const removed=await page.evaluate(()=>window.__ramDensityChecks.removed);
    assert.equal(removed,1);
    // Restore the original NGL Stage methods before changing structures.
    await page.evaluate(()=>{
      const stage=window.getNGLStage("NGL");
      if(stage && stage.__ramOriginalLoadFile)
        stage.loadFile=stage.__ramOriginalLoadFile;
      if(stage && stage.__ramOriginalRemoveComponent)
        stage.removeComponent=stage.__ramOriginalRemoveComponent;
    });

    // Replace the structure with a real NMR ensemble and run model analysis.
    await (await page.$("#structfile")).uploadFile(nmr);
    await new Promise(resolve=>setTimeout(resolve,1400));
    await page.click("#submit");
    await page.waitForFunction(()=>{
      const status=document.querySelector("#ram-current-structure");
      return status && status.textContent.includes("1D3Z");
    },{timeout:90000});
    // Shiny intentionally suspends outputs in hidden tabs, so activate
    // Summary before waiting for the lazily generated ensemble panel.
    await page.evaluate(()=>{
      const anchor=[...document.querySelectorAll(".nav-tabs a")]
        .find(x=>x.textContent.trim()==="Summary");
      if(!anchor)throw Error("Summary tab unavailable.");
      anchor.click();
    });
    await page.waitForFunction(()=>document.querySelector("#ram-ensemble-panel"),
                               {timeout:30000});
    await page.click("#ram-ensemble-panel summary");
    await page.click("#calculateEnsemble");
    try {
      await page.waitForFunction(()=>{
        // DataTables can create a separate header table when its tab is
        // resized; inspect rows across the output, not only its first table.
        const rows=document.querySelectorAll("#ensembleRows tbody tr");
        const summary=document.querySelector("#ensembleResultSummary");
        return !!(rows.length>2&&summary&&
          summary.textContent.includes("models analysed"));
      },{timeout:90000});
    }catch(error){
      const details=await page.evaluate(()=>({
        currentStructure:document.querySelector("#ram-current-structure")?.textContent,
        ensemblePanel:document.querySelector("#ram-ensemble-panel")?.textContent.slice(0,500),
        summary:document.querySelector("#ensembleResultSummary")?.textContent,
        tableRows:document.querySelectorAll("#ensembleRows table tbody tr").length,
        currentTab:document.querySelector(".nav-tabs li.active a")?.textContent,
        notifications:[...document.querySelectorAll(".shiny-notification")]
          .map(node=>node.textContent)
      }));
      console.error("Experimental verification ensemble diagnostics:",JSON.stringify(details));
      await page.screenshot({path:path.join(output,"verification-ensemble-failure.png"),
                             fullPage:true});
      throw error;
    }
    await page.screenshot({path:path.join(output,"verification-ensemble.png"),
                           fullPage:true});
    console.log("Experimental verification wwPDB import, NGL overlay and NMR ensemble passed.");
  }finally{await browser.close();}
})().catch(error=>{console.error(error);process.exitCode=1;});
