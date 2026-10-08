// Standalone lightweight in-app guide and responsive group layout regression.
// The application is already running at 127.0.0.1:8765 in ui-preview CI.
const puppeteer=require("puppeteer-core");
const fs=require("node:fs");
const path=require("node:path");
const assert=require("node:assert/strict");
(async()=>{
  const candidates=[
    process.env.CHROME_BIN,"/usr/bin/google-chrome",
    "/usr/bin/google-chrome-stable","/usr/bin/chromium",
    "/usr/bin/chromium-browser"].filter(Boolean);
  const executablePath=candidates.find(fs.existsSync);
  if(!executablePath)throw Error("No chromium installed.");
  const browser=await puppeteer.launch({
    executablePath,headless:true,
    args:["--no-sandbox","--disable-dev-shm-usage","--disable-gpu"]});
  const page=await browser.newPage();
  const errors=[];
  page.on("pageerror",e=>errors.push(String(e&&e.stack||e)));
  const dir=path.resolve("benchmarks/output/ui-preview");
  fs.mkdirSync(dir,{recursive:true});
  try{
    await page.setViewport({width:1440,height:850,deviceScaleFactor:1});
    await page.goto("http://127.0.0.1:8765",{
      waitUntil:"networkidle2",timeout:60000});
    await page.waitForSelector("#openGuide");
    await page.click("#openGuide");
    await page.waitForFunction(()=>!!document.querySelector(
      '.ram-main .nav-tabs li.active a[data-value="guide"]'),{timeout:12000});
    for(const id of ["ram-guide-inspect","ram-guide-prediction",
      "ram-guide-pair","ram-guide-groups","ram-guide-atlas",
      "ram-guide-export","ram-guide-interpretation"])
      assert.ok(await page.$("#"+id),"Missing guide content: "+id);
    const contentText=await page.$eval(".ram-guide",el=>el.textContent);
    for(const label of ["Rama8000","AlphaFold","Analyse groups",
      "Verify SIFTS mapping","30°","CSV"])
      assert.ok(contentText.includes(label),"Guide lacks: "+label);
    await page.screenshot({path:path.join(dir,"guide-desktop.png"),fullPage:true});
    await page.click("#guideGoGroups");
    await page.waitForFunction(()=>!!document.querySelector(
      '.ram-main .nav-tabs li.active a[data-value="compare"]'),{timeout:12000});
    await page.evaluate(()=>{
      document.getElementById("ram-group-comparison-panel").open=true;
    });
    await page.waitForSelector("#groupMinCoverage",{visible:true});
    assert.equal(await page.$("#downloadGroupComparison"),null,
      "Empty comparison should not offer an unusable CSV button.");
    async function checkLayout(width,height,name) {
      await page.setViewport({width,height,deviceScaleFactor:1});
      await new Promise(done=>setTimeout(done,250));
      const diagnostics=await page.evaluate(()=>{
        const root=document.querySelector(".ram-group-compare-controls");
        const cards=[...document.querySelectorAll(
          ".ram-group-compare-controls .ram-group-field")];
        const labelRects=[...document.querySelectorAll(
          ".ram-group-compare-controls label")].map(x=>{
          const r=x.getBoundingClientRect();return {left:r.left,right:r.right,
            top:r.top,bottom:r.bottom,text:x.innerText};});
        const fields=cards.map(x=>{
          const r=x.getBoundingClientRect();return {left:r.left,right:r.right,
            top:r.top,bottom:r.bottom,width:r.width};});
        const grid=root.getBoundingClientRect();
        const documentWidth=document.documentElement.clientWidth;
        const visible=fields.every(x=>x.width>=110&&x.left>=grid.left-2&&
          x.right<=grid.right+2);
        const overlaps=labelRects.some((x,i)=>labelRects.some((y,j)=>
          j>i && x.right>y.left+1 && x.left<y.right-1 &&
          x.bottom>y.top+1 && x.top<y.bottom-1));
        const uploadGrid=document.querySelector(".ram-group-upload-grid");
        const cardsRect=[...uploadGrid.children].map(x=>{
          const r=x.getBoundingClientRect();return {left:r.left,
            right:r.right,top:r.top,bottom:r.bottom};});
        const cardOverlaps=cardsRect.length>=2 &&
          cardsRect[0].left<cardsRect[1].right-1 &&
          cardsRect[0].right>cardsRect[1].left+1 &&
          cardsRect[0].top<cardsRect[1].bottom-1 &&
          cardsRect[0].bottom>cardsRect[1].top+1;
        return {fields,labels:labelRects,gridWidth:grid.width,visible,
          overlaps,cardOverlaps,
          docScrollWidth:document.documentElement.scrollWidth,
          documentWidth};
      });
      assert.ok(diagnostics.visible,
        "Group form fields escape container at "+width+": "+JSON.stringify(diagnostics));
      assert.equal(diagnostics.overlaps,false,
        "Group form labels overlap at "+width+": "+JSON.stringify(diagnostics));
      assert.equal(diagnostics.cardOverlaps,false,
        "Group upload cards overlap at "+width+": "+JSON.stringify(diagnostics));
      await page.screenshot({path:path.join(dir,name),fullPage:true});
    }
    await checkLayout(1440,900,"group-layout-desktop.png");
    await checkLayout(900,800,"group-layout-laptop.png");
    await checkLayout(390,800,"group-layout-mobile.png");

    await page.click("#groupGuide");
    await page.waitForFunction(()=>!!document.querySelector(
      '.ram-main .nav-tabs li.active a[data-value="guide"]'),{timeout:12000});
    await page.click("#guideGoAtlas");
    await page.waitForFunction(()=>!!document.querySelector(
      '.ram-main .nav-tabs li.active a[data-value="atlas"]'),{timeout:12000});
    await page.click("#atlasGuide");
    await page.waitForFunction(()=>!!document.querySelector(
      '.ram-main .nav-tabs li.active a[data-value="guide"]'),{timeout:12000});
    assert.equal(errors.length,0,"Browser JS errors: "+errors.join("\n"));
    console.log("Guide navigation and group layout passed at 1440, 900, 390 px.");
  }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
