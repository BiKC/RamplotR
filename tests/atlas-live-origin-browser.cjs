"use strict";
// This independent workflow loads the actual public Shinylive *origin* before
// issuing RCSB/PDBe fetches in Chromium. It does not require the newly built
// app to have been manually deployed: the PR's parser is injected from source.
// Mocked localhost browser tests cannot establish this external CORS property.
const assert=require("node:assert/strict");
const fs=require("node:fs");
const path=require("node:path");
const puppeteer=require("puppeteer-core");

(async()=>{
  const chrome=[process.env.CHROME_BIN,"/usr/bin/google-chrome",
    "/usr/bin/google-chrome-stable","/usr/bin/chromium"].find(
      v=>v && fs.existsSync(v));
  if(!chrome) throw Error("Chromium executable is required.");
  const browser=await puppeteer.launch({executablePath:chrome,headless:true,
    args:["--no-sandbox","--disable-dev-shm-usage","--disable-gpu"]});
  const output="benchmarks/output/atlas-live-browser-origin";
  fs.mkdirSync(output,{recursive:true});
  try{
    const page=await browser.newPage();
    page.setDefaultTimeout(45000);
    const errors=[];
    page.on("pageerror",e=>errors.push(String(e.message||e).slice(0,250)));
    const response=await page.goto("https://bikc.be/RamplotR/",{
      waitUntil:"domcontentloaded",timeout:60000
    });
    const origin=await page.evaluate(()=>window.location.origin);
    assert.equal(origin,"https://bikc.be",
      "Deployed browser origin is not https://bikc.be.");
    assert.ok(response && response.status()>=200 && response.status()<400,
      "Public Shinylive page returned an unsuccessful HTTP response.");
    for(const file of ["atlas-discovery.js","atlas-sifts-exact.js",
                       "atlas-connectivity.js"])
      await page.addScriptTag({path:path.join("shinyRam","www",file)});
    const results=await page.evaluate(async()=>
      window.RamplotRAtlasConnectivity.runChecks({timeoutMs:22000}));
    fs.writeFileSync(path.join(output,"origin-checks.json"),
      JSON.stringify({page_url:page.url(),http_status:response.status(),
        ...results,page_errors:errors},null,2)+"\n");
    console.log("Live deployed-origin RCSB/PDBe browser diagnostic:");
    console.log(JSON.stringify(results,null,2));
    assert.equal(results.tested_origin,origin);
    assert.equal(results.passed,3,
      "The public-site browser origin cannot access all three endpoints.");
  }finally{await browser.close();}
})().catch(error=>{console.error(error.stack||error);process.exitCode=1;});
