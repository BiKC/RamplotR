// Test PAE heatmap linkage without a browser or an R/Shiny session.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const vm = require("node:vm");

(async () => {
  const handlers = {};
  const events = {};
  const clicks = {};
  const inputs = [];
  let draw = null;
  let restyle = null;
  const note = {textContent: ""};
  const panel = {open: true};
  const node = {
    id: "ram-pae-plot",
    parentElement: {clientWidth: 520},
    data: null,
    on(name, callback) { clicks[name] = callback; }
  };
  const document = {
    getElementById(id) {
      return {"ram-pae-plot":node, "ram-confidence-panel":panel,
              "ram-pae-note":note}[id] || null;
    },
    addEventListener(name, callback) { events[name] = callback; }
  };
  const window = {
    requestAnimationFrame(callback) { callback(); },
    addEventListener() {},
    Shiny: {
      addCustomMessageHandler(name, callback) {
        assert.equal(callback.length, 1, name + " must accept a Shiny message");
        handlers[name] = callback;
      },
      setInputValue(name, value) { inputs.push([name, value]); }
    },
    Plotly: {
      react(target, traces, layout, config) {
        draw = {target, traces, layout, config};
        node.data = traces;
        return Promise.resolve();
      },
      restyle(target, changes, traces) {
        restyle = {target, changes, traces};
      }
    }
  };
  vm.runInNewContext(
    fs.readFileSync("shinyRam/www/prediction.js", "utf8"),
    {window, document});
  assert.equal(typeof handlers["ram-confidence"], "function");
  assert.equal(typeof handlers["ram-confidence-selected"], "function");
  const residues = [1, 2, 3].map(resi => ({
    chain:"A",resi,insertion_code:""
  }));
  handlers["ram-confidence"]({
    z: [[0,4,8],[2,0,5],[7,3,0]],
    labels:["A:1","A:2","A:3"],
    residues,total_tokens:3,displayed_tokens:3,downsampled:false
  });
  await new Promise(resolve => setImmediate(resolve));
  assert.equal(draw.target, node);
  assert.equal(draw.traces[0].z.length, 3);
  assert.equal(draw.layout.yaxis.scaleanchor, "x");
  assert.equal(typeof clicks.plotly_click, "function");
  clicks.plotly_click({points:[{x:2,y:1}]});
  assert.equal(inputs.length, 1);
  assert.equal(inputs[0][0], "ramPaePick");
  assert.equal(inputs[0][1], residues[1]);

  handlers["ram-confidence-selected"](residues[2]);
  assert.equal(restyle.target,node);
  assert.equal(restyle.changes.x[0][0],3);
  assert.equal(restyle.changes.y[0][0],3);
  handlers["ram-confidence"]({clear:true});
  await new Promise(resolve => setImmediate(resolve));
  assert.equal(note.textContent.includes("directional"), true);
  console.log("Prediction confidence plot linkage tests passed");
})().catch(error => { console.error(error); process.exitCode = 1; });
