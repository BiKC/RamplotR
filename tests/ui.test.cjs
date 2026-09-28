// Node-only behaviour checks. No browser, network or Shiny server required.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const vm = require("node:vm");

const handlers = {};
const listeners = {};
let selected = "pdb";
let rendered = null;
const events = {};
const inputs = [];
const restyles = [];
let nglClick = null;
let nglReady = false;
const stickSelections = [];
const cameraMoves = [];
const rockEvents = [];
const spinEvents = [];
const stage = {
  signals: { clicked: {
    add(fn) { nglClick = fn; },
    remove(fn) { if (nglClick === fn) nglClick = null; }
  }},
  getRepresentationsByName(name) {
    assert.equal(name, "ram-highlight");
    return { list: [{}], setSelection(sele) { stickSelections.push(sele); } };
  },
  setRock(value) { rockEvents.push(value); },
  setSpin(value) { spinEvents.push(value); },
  autoView(duration) { cameraMoves.push({ overview: true, duration }); }
};
const structure = {
  structure: {},
  autoView(sele, duration) { cameraMoves.push({ sele, duration }); }
};

function fakeClassList() {
  const names = new Set();
  return {
    toggle(name, condition) {
      if (condition) names.add(name);
      else names.delete(name);
    },
    contains(name) { return names.has(name); }
  };
}

const pdbWrap = { classList: fakeClassList() };
const uploadWrap = { classList: fakeClassList() };
const plot = {
  clientWidth: 650,
  style: {},
  classList: { contains: () => true },
  on(name, fn) { events[name] = fn; }
};
const placeholderText = { textContent: "" };
const empty = {
  hidden: false,
  querySelector: () => placeholderText
};
const badge = { textContent: "" };
const elements = {
  plotly: plot,
  "plot-empty": empty,
  "ram-current-structure": badge,
  "ram-pdb-wrap": pdbWrap,
  "ram-upload-wrap": uploadWrap
};
const document = {
  readyState: "complete",
  getElementById(id) { return elements[id] || null; },
  querySelector(selector) {
    return selector === 'input[name="inputSource"]:checked'
      ? { value: selected } : null;
  },
  addEventListener(name, callback) { listeners[name] = callback; }
};
const window = {
  Shiny: {
    addCustomMessageHandler(name, callback) { handlers[name] = callback; },
    setInputValue(name, value) { inputs.push({name, value}); }
  },
  Plotly: {
    react(node, traces, layout, config) {
      rendered = { node, traces, layout, config };
      return Promise.resolve();
    },
    restyle(node, change, traces) { restyles.push({node, change, traces}); },
    Plots: { resize() {} }
  }
};
window.getNGLStage = () => nglReady ? stage : null;
window.getNGLStructure = () => nglReady ? [structure] : null;
const source = fs.readFileSync("shinyRam/www/custom.js", "utf8");
vm.runInNewContext(source, { window, document, Event });

assert.equal(typeof handlers.process, "function",
             "Shiny plot message handler must be registered");
assert.equal(uploadWrap.classList.contains("is-hidden"), true,
             "PDB input should be visible initially");
selected = "upload";
listeners.change({ target: { name: "inputSource" } });
assert.equal(pdbWrap.classList.contains("is-hidden"), true);
assert.equal(uploadWrap.classList.contains("is-hidden"), false);

handlers.process({
  name: "1CRN",
  df: {
    chain: ["A", "A"],
    resn: ["ALA", "<GLY>"],
    resi: [1, 2],
    insertion_code: ["", "A"],
    phi: [null, 91.2],
    psi: [null, -37.1]
  },
  matrix: {
    x: [-180, 0, 180],
    y: [-180, 0, 180],
    z: [[1, 2, 3], [3, 4, 5], [5, 6, 7]]
  },
  limits: [5, 3, 1],
  backgroundColors: ["#F1EEF6", "#BDC9E1", "#74A9CF", "#0570B0"],
  chainColors: ["#116e70"]
});
assert.equal(rendered.node, plot, "The expected plot container must be used");
assert.equal(rendered.traces.length, 6, "Four contours, one chain and one selection trace");
assert.equal(rendered.traces[4].x.length, 1,
             "A missing torsion must not turn into a false (0,0) point");
assert.equal(rendered.traces[4].x[0], 91.2);
assert.ok(rendered.traces[4].text[0].includes("&lt;GLY&gt;"),
          "Hover labels must escape structure-supplied strings");
assert.equal(rendered.layout.yaxis.scaleanchor, "x",
             "Angular axes must retain the same scale");
assert.equal(rendered.config.responsive, true);
assert.equal(empty.hidden, true);
assert.ok(badge.textContent.includes("1 plotted residues"));
// Plotly.react promises resolve on the next microtask; point-click handlers
// are registered after rendering, not before.
setImmediate(() => {
  try {
  assert.equal(typeof events.plotly_click, "function");
  assert.equal(rendered.traces[4].customdata[0][2], "A",
               "Insertions must stay attached to the plotted residue.");
  // NGL may finish loading after the 2D plot was clicked.
  handlers["ram-selection"]({ chain: "A", resi: 2, insertion_code: "A" });
  assert.equal(restyles.length, 1, "Linked selection must highlight the plot.");
  assert.equal(restyles[0].change.x[0][0], 91.2);
  assert.equal(stickSelections.length, 0,
               "Do not manipulate NGL before it is ready.");
  nglReady = true;
  handlers["ram-bind-ngl"]({});
  assert.equal(stickSelections.at(-1), "2^A:A",
               "A selected residue must be shown in orange sticks.");
  assert.deepEqual(cameraMoves.at(-1), { sele: "2^A:A", duration: 650 },
                   "NGL should animate camera focus to the selected residue.");
  assert.deepEqual(rockEvents, [false], "Stop rocking during residue zoom.");
  assert.deepEqual(spinEvents, [false], "Stop spinning during residue zoom.");
  handlers["ram-selection"]({ chain: "A", resi: 2, insertion_code: "A" });
  assert.equal(cameraMoves.length, 1,
               "Repeated plot redraws must not refocus the same residue.");
  handlers["ram-selection"]({ chain: "A", resi: 1, insertion_code: "" });
  assert.deepEqual(cameraMoves.at(-1), { sele: "1:A", duration: 650 },
                   "Selecting a different residue must move the camera.");
  handlers["ram-selection"]({ chain: "A", resi: 2, insertion_code: "A" });
  events.plotly_click({
    points: [{ customdata: ["A", 2, "A", "GLY"] }]
  });
  assert.equal(inputs.at(-1).name, "ramPlotPick");
  assert.equal(inputs.at(-1).value.insertion_code, "A");
  handlers["ram-bind-ngl"]({});
  assert.equal(typeof nglClick, "function",
               "The NGL stage must listen for atom picking.");
  nglClick({ atom: { chainname: "A", resno: 2, inscode: "A" } });
  assert.equal(inputs.at(-1).name, "ramNglPick",
               "An NGL click must publish the same residue identity.");
  assert.equal(inputs.at(-1).value.resi, 2);
  handlers["ram-selection"]({clear: true});
  assert.equal(restyles.at(-1).change.x[0].length, 0,
               "Clear selection must remove the overlay.");
  assert.equal(stickSelections.at(-1), "none",
               "Clear selection must hide residue sticks.");
  assert.deepEqual(cameraMoves.at(-1), { overview: true, duration: 450 },
                   "Clearing selection must restore the whole structure view.");
  console.log("RamplotR plot, source switching and three-way picking tests passed.");
  } catch (e) { console.error(e); process.exitCode = 1; }
});
