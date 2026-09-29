// No network required: verify that Plotly downloads only after it is needed.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const vm = require("node:vm");

const source = fs.readFileSync("shinyRam/www/plotly-loader.js", "utf8");

function setup() {
  const scripts = [];
  const events = {};
  const document = {
    createElement(tag) {
      assert.equal(tag, "script");
      return { parentNode: { removeChild() {} } };
    },
    head: { appendChild(script) { scripts.push(script); } },
    addEventListener(name, callback) { events[name] = callback; }
  };
  const window = {};
  vm.runInNewContext(source, { window, document, Promise, Error });
  return { scripts, events, window };
}

(async function () {
  const app = setup();
  assert.equal(app.scripts.length, 0, "No Plotly request on initial load");
  assert.equal(typeof app.window.ramLoadPlotly, "function");
  const first = app.window.ramLoadPlotly();
  const second = app.window.ramLoadPlotly();
  assert.equal(first, second, "Concurrent plot requests share one download");
  assert.equal(app.scripts.length, 1);
  assert.match(app.scripts[0].src, /plotly-2\.14\.0\.min\.js$/);
  app.window.Plotly = { react() {} };
  app.scripts[0].onload();
  assert.equal(await first, app.window.Plotly);
  assert.equal(await app.window.ramLoadPlotly(), app.window.Plotly);
  assert.equal(app.scripts.length, 1, "Cached Plotly needs no new request");

  const early = setup();
  early.events.click({target: {closest() {return {id: "submit"};}}});
  assert.equal(early.scripts.length, 1, "Analyse starts Plotly download");
  early.window.Plotly = {};
  early.scripts[0].onload();
  await early.window.ramLoadPlotly();

  const failed = setup();
  const request = failed.window.ramLoadPlotly();
  failed.scripts[0].onerror();
  await assert.rejects(request, /Could not download Plotly/);
  const retry = failed.window.ramLoadPlotly();
  assert.equal(failed.scripts.length, 2, "Failed downloads can be retried");
  failed.window.Plotly = {};
  failed.scripts[1].onload();
  await retry;
  console.log("Deferred Plotly loader tests passed.");
})().catch(error => {
  console.error(error);
  process.exitCode = 1;
});
