// Node-only behaviour checks. No browser, network or Shiny server required.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const vm = require("node:vm");

const handlers = {};
const listeners = {};
let selected = "pdb";
let rendered = null;

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
  classList: { contains: () => true }
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
    addCustomMessageHandler(name, callback) { handlers[name] = callback; }
  },
  Plotly: {
    react(node, traces, layout, config) {
      rendered = { node, traces, layout, config };
      return Promise.resolve();
    },
    Plots: { resize() {} }
  }
};
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
assert.equal(rendered.traces.length, 5, "Four contours plus one chain trace");
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
console.log("RamplotR JS source controls and plot rendering checks passed.");
