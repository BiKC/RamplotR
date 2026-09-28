// Browser-independent NGL volume handler contract; no actual experimental data.
const assert=require("node:assert/strict");
const vm=require("node:vm"),fs=require("node:fs");
(async function(){
  const listeners={};const handlers={};
  const nodes={
    "ram-density-file":{files:[{name:"local-map.mrc",size:2048}]},
    "ram-density-level":{value:"2"},
    "ram-density-value":{textContent:""},
    "ram-density-status":{textContent:""}
  };
  const events={
    load:{target:{closest:selector=>selector==="#ram-density-load"}},
    clear:{target:{closest:selector=>selector==="#ram-density-clear"}}
  };
  const representation={changes:[],setParameters(p){this.changes.push(p);}};
  const component={reprList:[representation],
                   addRepresentation(type,params) {
                     assert.equal(type,"surface");
                     assert.equal(params.isolevelType,"sigma");
                     return representation;
                   }};
  const removed=[];
  const stage={loadFile:async(file,options)=>{
      assert.equal(file.name,"local-map.mrc");
      assert.equal(options.ext,"ccp4");
      return component;
    },removeComponent(value){removed.push(value);}};
  const sandbox={
    document:{
      getElementById:id=>nodes[id]||null,
      addEventListener:(event,handler)=>listeners[event]=handler
    },
    window:{
      getNGLStage:()=>stage,
      Shiny:{addCustomMessageHandler:(name,handler)=>{
        assert.equal(handler.length,1,"Shiny handlers require one message argument");
        handlers[name]=handler;
      }}
    },console,Promise,Number
  };
  vm.runInNewContext(fs.readFileSync("shinyRam/www/density.js","utf8"),
                     sandbox);
  assert.equal(typeof handlers["ram-clear-density"],"function");
  listeners.click(events.load);
  await new Promise(resolve=>setImmediate(resolve));
  assert.match(nodes["ram-density-status"].textContent,/2.00σ/);
  nodes["ram-density-level"].value="3.25";
  listeners.change({target:{id:"ram-density-level"}});
  assert.equal(representation.changes.length,1);
  assert.equal(representation.changes[0].isolevel,3.25);
  assert.match(nodes["ram-density-value"].textContent,/3.25σ/);
  handlers["ram-clear-density"]({});
  assert.equal(removed.length,1);
  assert.match(nodes["ram-density-status"].textContent,/No map/);
  nodes["ram-density-file"].files=[{name:"huge.ccp4",size:70*1024*1024}];
  listeners.click(events.load);
  await new Promise(resolve=>setImmediate(resolve));
  assert.match(nodes["ram-density-status"].textContent,/64 MB/);
  assert.equal(removed.length,1);
  console.log("Local NGL map controls and safety checks passed.");
})().catch(error=>{console.error(error);process.exitCode=1;});
