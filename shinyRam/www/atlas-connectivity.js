(function (root) {
  "use strict";

  // Run only on explicit user request. All three probes execute directly in
  // the visitor's browser (including Shinylive/webR), not from an R process,
  // Actions runner, proxy or different origin.
  const SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query";
  const META_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/4AKE/1";
  const SIFTS_URL = "https://www.ebi.ac.uk/pdbe/static/entry/4ake_updated.cif";
  const TEST_PDB = "4AKE", TEST_ACCESSION = "P69441";
  const MAX_MESSAGE = 240;
  let activeRequest = 0;

  function notify(value) {
    if(root.Shiny && typeof root.Shiny.setInputValue === "function")
      root.Shiny.setInputValue("ramAtlasConnectivity",value,{priority:"event"});
  }
  function diagnostic(error) {
    const text = String(error && error.message || error || "Unknown error")
      .replace(/\s+/g," ").slice(0,MAX_MESSAGE);
    const timedOut = error && error.name === "AbortError";
    const network = error && error.name === "TypeError";
    return {
      state:"failed",
      reason:timedOut?"timeout":network?"network_or_cors":"invalid_response",
      detail: timedOut ? "Request timed out." :
        network ? "Browser could not fetch this endpoint. Check network, CORS policy, firewall or privacy extensions." :
        text
    };
  }
  async function request(fetcher,url,options={},timeoutMs=15000) {
    const controller = typeof AbortController === "function"
      ? new AbortController() : null;
    let timer=null;
    if(controller) timer=setTimeout(()=>controller.abort(),timeoutMs);
    try {
      const response=await fetcher(url,{
        ...options,credentials:"omit",cache:"no-store",
        signal:controller ? controller.signal : undefined
      });
      if(!response || typeof response.status !== "number")
        throw new Error("Not a valid browser HTTP response.");
      if(!response.ok) {
        const error=new Error("HTTP "+response.status+
          ". The archive returned a non-success response.");
        error.httpStatus=response.status;
        throw error;
      }
      return response;
    } finally {
      if(timer!==null) clearTimeout(timer);
    }
  }
  function result(id,name,url,start,extra) {
    const elapsed=Math.max(0,Math.round(Date.now()-start));
    return {id,name,url,elapsed_ms:elapsed,...extra};
  }
  const descriptors=[
    {id:"rcsb_search",name:"RCSB experimental entity search",url:SEARCH_URL},
    {id:"rcsb_metadata",name:"RCSB polymer entity metadata",url:META_URL},
    {id:"pdbe_sifts",name:"PDBe updated mmCIF / exact SIFTS",url:SIFTS_URL}
  ];

  async function check(item,fetcher,timeoutMs) {
    const start=Date.now();
    try {
      let response,summary;
      if(item.id==="rcsb_search") {
        const api=root.RamplotRAtlasDiscovery;
        if(!api || typeof api.searchRequest!=="function")
          throw new Error("RamplotR Atlas search parser is unavailable.");
        response=await request(fetcher,item.url,{
          method:"POST",
          headers:{"Content-Type":"application/json"},
          body:JSON.stringify(api.searchRequest(TEST_ACCESSION,1,0))
        },timeoutMs);
        const data=await response.json();
        if(!data || !Array.isArray(data.result_set) ||
           !Number.isFinite(Number(data.total_count)))
          throw new Error("RCSB response lacks expected search results.");
        summary=String(data.total_count)+" experimental polymer-entity hits";
      } else if(item.id==="rcsb_metadata") {
        response=await request(fetcher,item.url,{},timeoutMs);
        const data=await response.json();
        if(!data || !data.rcsb_polymer_entity_container_identifiers)
          throw new Error("RCSB returned an unexpected polymer-entity payload.");
        summary="4AKE entity 1 metadata parsed";
      } else if(item.id==="pdbe_sifts") {
        const api=root.RamplotRExactSifts;
        if(!api || typeof api.exactRows!=="function" ||
           typeof api.extractMappedCA!=="function")
          throw new Error("RamplotR SIFTS parser is unavailable.");
        response=await request(fetcher,item.url,{},timeoutMs);
        const text=await response.text();
        if(text.length>45*1024*1024)
          throw new Error("Updated mmCIF is larger than the 45 MB safety limit.");
        const exact=api.exactRows(text,TEST_PDB,"1",TEST_ACCESSION);
        const ca=api.extractMappedCA(text,exact);
        const observed=exact.rows.filter(x=>x.observed).length;
        if(observed<1 || ca.length<1)
          throw new Error("No observed SIFTS-mapped C-alpha residues were parsed.");
        summary=observed+" observed mapped residue rows; "+
          ca.length+" first-model C-alpha atoms";
      } else {
        throw new Error("Unsupported Atlas connectivity probe.");
      }
      return result(item.id,item.name,item.url,start,
        {state:"ok",http_status:response.status,summary});
    } catch(error) {
      return result(item.id,item.name,item.url,start,{
        ...diagnostic(error),
        ...(Number.isInteger(error && error.httpStatus)
            ? {http_status:error.httpStatus}:{})
      });
    }
  }
  async function runChecks(options={}) {
    const fetcher=options.fetcher || (typeof root.fetch==="function"
      ? root.fetch.bind(root) : null);
    const limit=Number(options.timeoutMs);
    const timeoutMs=Number.isFinite(limit)&&limit>=500&&limit<=60000
      ? limit : 15000;
    const rows=[];
    if(!fetcher) throw new Error("Browser fetch API is unavailable.");
    for(const item of descriptors) {
      rows.push(await check(item,fetcher,timeoutMs));
      if(typeof options.onProgress==="function")
        options.onProgress(rows.slice());
    }
    const passed=rows.filter(x=>x.state==="ok").length;
    return {
      state:passed===rows.length?"ok":"partial_failure",
      tested_origin:root.location && root.location.origin ||
        "unknown browser origin",
      requested_at:new Date().toISOString(),
      passed,total:rows.length,checks:rows
    };
  }
  async function runForShiny(payload) {
    const id=String(payload && payload.request_id || "");
    const seq=++activeRequest;
    notify({request_id:id,state:"running",checks:[],passed:0,total:3});
    try {
      const value=await runChecks({onProgress:rows=>{
        if(seq!==activeRequest) return;
        notify({request_id:id,state:"running",checks:rows,
          passed:rows.filter(x=>x.state==="ok").length,total:3});
      }});
      if(seq===activeRequest) notify({request_id:id,...value});
    } catch(error) {
      if(seq===activeRequest) notify({request_id:id,state:"error",
        message:diagnostic(error).detail,checks:[],passed:0,total:3});
    }
  }
  function setup() {
    if(!root.Shiny ||
       typeof root.Shiny.addCustomMessageHandler!=="function")return false;
    root.Shiny.addCustomMessageHandler("ram-atlas-connectivity",runForShiny);
    return true;
  }
  if(!setup() && root.document)
    root.document.addEventListener("shiny:connected",setup,{once:true});
  const api={descriptors,diagnostic,check,runChecks};
  root.RamplotRAtlasConnectivity=api;
  if(typeof module!=="undefined" && module.exports)module.exports=api;
})(typeof window!=="undefined" ? window : globalThis);
