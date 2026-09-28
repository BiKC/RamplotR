// Run on CI with a real Shiny server and Chromium. Archive screenshots to
// inspect desktop and mobile layouts without deploying the branch.
const puppeteer = require("puppeteer-core");
const path = require("node:path");
const fs = require("node:fs");
const assert = require("node:assert/strict");

(async function () {
  const candidates = [
    process.env.CHROME_BIN,
    "/usr/bin/google-chrome",
    "/usr/bin/google-chrome-stable",
    "/usr/bin/chromium",
    "/usr/bin/chromium-browser"
  ].filter(Boolean);
  const executablePath = candidates.find(fs.existsSync);
  if (!executablePath) throw Error("Chromium or Chrome was not found.");
  const browser = await puppeteer.launch({
    executablePath, headless: true,
    args: ["--no-sandbox", "--disable-gpu", "--disable-dev-shm-usage"]
  });
  const errors = [];
  try {
    const page = await browser.newPage();
    page.on("pageerror", error => errors.push(String(error && (error.stack || error.message || error))));
    await page.setViewport({ width: 1440, height: 940, deviceScaleFactor: 1 });
    await page.goto("http://127.0.0.1:8765", {
      waitUntil: "networkidle2", timeout: 60000
    });
    await page.waitForSelector(".ram-workspace");
    assert.ok(await page.$("#ram-pdb-wrap"));
    assert.ok(await page.$("#ram-upload-wrap"));
    assert.ok(await page.$("#NGL"));
    await page.screenshot({
      path: "benchmarks/output/ui-preview/desktop-empty.png", fullPage: true
    });

    // Upload a local real protein to avoid depending on an RCSB fetch inside
    // the Shiny session. It exercises both the source selector and plotting.
    await page.evaluate(() => {
      document.querySelector('input[name="inputSource"][value="upload"]').click();
    });
    await page.waitForFunction(() =>
      !document.getElementById("ram-upload-wrap").classList.contains("is-hidden"));
    const input = await page.$("#structfile");
    await input.uploadFile(path.resolve("benchmarks/output/ui-preview/1CRN.pdb"));
    // Shiny's upload widget does not expose a stable progress-complete
    // attribute across versions. Confirm the file was selected, then give the
    // small local upload a moment to finish before submitting.
    await page.waitForFunction(() => {
      const field = document.querySelector("#structfile");
      return field && field.files && field.files.length === 1;
    }, { timeout: 15000 });
    await new Promise(resolve => setTimeout(resolve, 1500));
    await page.click("#submit");
    await page.waitForFunction(() =>
      document.getElementById("plotly").classList.contains("js-plotly-plot"),
      { timeout: 90000 }
    );
    const plot = await page.evaluate(() => {
      const p = document.getElementById("plotly");
      return {
        traces: p.data.length,
        yAnchor: p._fullLayout.yaxis.scaleanchor,
        xRange: p._fullLayout.xaxis.range,
        yRange: p._fullLayout.yaxis.range,
        axisPixels: {
          x: p._fullLayout.xaxis._length,
          y: p._fullLayout.yaxis._length
        },
        renderedHeight: p.getBoundingClientRect().height,
        status: document.getElementById("ram-current-structure").textContent,
        width: p.getBoundingClientRect().width
      };
    });
    assert.ok(plot.traces >= 5, "Expected contours and at least one chain.");
    assert.equal(plot.yAnchor, "x", "Axes must be equally scaled.");
    assert.ok(Math.abs(plot.xRange[0] + 180) < 1 &&
              Math.abs(plot.xRange[1] - 180) < 1 &&
              Math.abs(plot.yRange[0] + 180) < 1 &&
              Math.abs(plot.yRange[1] - 180) < 1,
              "The plot must keep both -180° to 180° angular ranges.");
    assert.ok(plot.status.includes("plotted residues"), "No plotted residues");
    assert.ok(plot.axisPixels.x >= 0.48 * plot.width &&
              plot.axisPixels.y >= 0.48 * plot.width,
              "Angular axes must use at least half of the plot panel width.");
    // NGL loads and paints asynchronously after the Plotly response.
    // Give the viewer a moment to render before taking the desktop preview.
    await new Promise(resolve => setTimeout(resolve, 1600));
    await page.screenshot({
      path: "benchmarks/output/ui-preview/desktop-loaded.png", fullPage: true
    });


    // Changing presentation controls must update the displayed plot without
    // another click on Analyze, and without loading another structure.
    await page.evaluate(() => {
      // selectInput is Selectize-backed; changing its hidden native select
      // bypasses the Shiny binding. Use the same API as a real user choice.
      const el = document.querySelector("#colorscheme");
      if (el.selectize) el.selectize.setValue("PDBSum");
      else {
        el.value = "PDBSum";
        el.dispatchEvent(new Event("change", { bubbles: true }));
      }
    });
    try {
      await page.waitForFunction(() => {
        const p = document.getElementById("plotly");
        return p && p.data && p.data[0].fillcolor === "#F3F300";
      }, { timeout: 18000 });
    } catch (e) {
      console.error("Palette diagnostics:", JSON.stringify(
        await page.evaluate(() => ({
          selected: document.querySelector("#colorscheme").value,
          pickerColors: ["bg1", "bg2", "bg3", "bg4"].map(id => {
            const el = document.getElementById(id);
            return el ? el.value : null;
          }),
          actualPlotColor: document.getElementById("plotly").data[0].fillcolor,
          currentStructure: document.getElementById("ram-current-structure").textContent
        }))
      ));
      await page.screenshot({
        path: "benchmarks/output/ui-preview/palette-debug.png", fullPage: true
      });
      throw e;
    }
    const fullCount = await page.evaluate(() => {
      const p = document.getElementById("plotly");
      return p.data.filter(trace => trace.customdata)
        .reduce((sum, trace) => sum + trace.x.length, 0);
    });
    assert.ok(fullCount > 5, "The initial test structure should contain many points.");
    await page.evaluate(() => window.Shiny.setInputValue("AA", ["GLY"], {
      priority: "event"
    }));
    await page.waitForFunction(expected => {
      const p = document.getElementById("plotly");
      const count = p && p.data && p.data.filter(t => t.customdata)
        .reduce((n, t) => n + t.x.length, 0);
      return count > 0 && count < expected;
    }, { timeout: 18000 }, fullCount);
    await page.evaluate(() => {
      const values = Array.from(document.querySelector("#AA").options)
        .map(option => option.value);
      window.Shiny.setInputValue("AA", values, { priority: "event" });
    });
    await page.waitForFunction(expected => {
      const p = document.getElementById("plotly");
      const count = p && p.data && p.data.filter(t => t.customdata)
        .reduce((n, t) => n + t.x.length, 0);
      return count === expected;
    }, { timeout: 18000 }, fullCount);

    // Instrument the actual NGL widget: a plot click must turn on the
    // persistent ball-and-stick representation AND animate camera focus.
    await page.waitForFunction(() => window.getNGLStage &&
      window.getNGLStructure && window.getNGLStructure("NGL") &&
      window.getNGLStructure("NGL").length &&
      window.getNGLStage("NGL").getRepresentationsByName("ram-highlight")
        .list.length > 0, { timeout: 25000 });
    await page.evaluate(() => {
      const stage = window.getNGLStage("NGL");
      const component = window.getNGLStructure("NGL")[0];
      window.__ramZoomCalls = [];
      window.__ramStickSelections = [];
      window.__ramOverviewCalls = [];
      const autoView = component.autoView.bind(component);
      component.autoView = (sele, duration) => {
        window.__ramZoomCalls.push({ sele, duration });
        return autoView(sele, duration);
      };
      const overview = stage.autoView.bind(stage);
      stage.autoView = duration => {
        window.__ramOverviewCalls.push(duration);
        return overview(duration);
      };
      const getRepresentations = stage.getRepresentationsByName.bind(stage);
      stage.getRepresentationsByName = name => {
        const group = getRepresentations(name);
        if (name === "ram-highlight" && group && group.setSelection) {
          const setSelection = group.setSelection.bind(group);
          group.setSelection = sele => {
            window.__ramStickSelections.push(sele);
            return setSelection(sele);
          };
        }
        return group;
      };
    });
    const firstPoint = await page.evaluate(() => {
      const trace = document.getElementById("plotly").data.find(t =>
        t.customdata && t.customdata.length);
      return trace.customdata[0];
    });
    await page.evaluate(detail => {
      document.getElementById("plotly").emit("plotly_click", {
        points: [{ customdata: detail }]
      });
    }, firstPoint);
    const firstNglSelector = firstPoint[1] +
      (firstPoint[2] ? "^" + firstPoint[2] : "") +
      (firstPoint[0] ? ":" + firstPoint[0] : "");
    await page.waitForFunction(expected => {
      const p = document.getElementById("plotly");
      return document.querySelector("#selectedResidueInfo strong") &&
        p.data[p.data.length - 1].x.length === 1 &&
        window.__ramStickSelections.includes(expected) &&
        window.__ramZoomCalls.some(call =>
          call.sele === expected && call.duration === 650);
    }, { timeout: 20000 }, firstNglSelector);
    await page.screenshot({
      path: "benchmarks/output/ui-preview/desktop-residue-zoom.png",
      fullPage: true
    });

    // The residue table should also select the same entry, and clicking
    // another row should move the plot highlight back to that residue.
    await page.click('.nav-tabs a[data-value="residues"]');
    await page.waitForSelector("#regions table tbody tr", { timeout: 18000 });
    await page.waitForFunction(() =>
      document.querySelector("#regions tbody tr.selected"), { timeout: 18000 });
    const before = await page.$eval("#selectedResidueInfo strong",
                                   e => e.textContent);
    // A single page.$ returns only the first row; find and actually click a
    // different visible residue so the table -> NGL path is exercised.
    // DT listens to actual pointer clicks on table cells. Calling tr.click()
    // dispatches an artificial row event that DT does not report to Shiny.
    const targetHandle = await page.evaluateHandle(first => {
      const rows = Array.from(document.querySelectorAll("#regions table tbody tr"));
      const row = rows.find(candidate => {
        const cells = candidate.querySelectorAll("td");
        return cells.length >= 3 &&
          (cells[0].textContent.trim() !== String(first[0]) ||
           Number(cells[1].textContent.trim()) !== Number(first[1]) ||
           cells[2].textContent.trim() !== String(first[2] || ""));
      });
      return row ? row.querySelector("td") : null;
    }, firstPoint);
    const targetCell = targetHandle.asElement();
    assert.ok(targetCell, "Need at least two different visible residues.");
    const other = await page.evaluate(cell => {
      const cells = cell.closest("tr").querySelectorAll("td");
      const chain = cells[0].textContent.trim();
      const resi = cells[1].textContent.trim();
      const insertion = cells[2].textContent.trim();
      return { sele: resi + (insertion ? "^" + insertion : "") +
                      (chain ? ":" + chain : "") };
    }, targetCell);
    await targetCell.click();
    await targetHandle.dispose();
    try {
      await page.waitForFunction(({ before, sele }) =>
        document.querySelector("#selectedResidueInfo strong") &&
        document.querySelector("#selectedResidueInfo strong").textContent !== before &&
        window.__ramStickSelections.includes(sele) &&
        window.__ramZoomCalls.some(call => call.sele === sele),
        { timeout: 20000 }, { before, sele: other.sele });
    } catch (error) {
      const diagnostic = await page.evaluate(({before, expected}) => ({
        before: before,
        expected: expected,
        selectedInfo: document.querySelector("#selectedResidueInfo")?.innerText,
        chosenRows: Array.from(document.querySelectorAll("#regions tbody tr.selected"))
          .map(row => row.innerText),
        dtLastClicked: window.Shiny?.shinyapp?.$inputValues?.["regions_row_last_clicked"],
        dtRowsSelected: window.Shiny?.shinyapp?.$inputValues?.["regions_rows_selected"],
        zoomCalls: window.__ramZoomCalls,
        stickSelections: window.__ramStickSelections,
        selectionOutput: document.getElementById("ram-current-structure")?.innerText
      }), {before, expected:other.sele});
      console.error("Table-to-NGL diagnostics:", JSON.stringify(diagnostic));
      console.error("Browser page errors:", errors);
      await page.screenshot({
        path: "benchmarks/output/ui-preview/table-selection-debug.png",
        fullPage: true
      });
      throw error;
    }
    await page.click('.nav-tabs a[data-value="plot"]');
    await page.click("#clearResidue");
    try {
      await page.waitForFunction(() => {
        const p = document.getElementById("plotly");
        return !document.querySelector("#selectedResidueInfo strong") &&
          p.data[p.data.length - 1].x.length === 0 &&
          window.__ramStickSelections.includes("none") &&
          window.__ramOverviewCalls.includes(450);
      }, { timeout: 12000 });
    } catch (e) {
      const state = await page.evaluate(() => {
        const p = document.getElementById("plotly");
        return {
          selectedInfo: document.querySelector("#selectedResidueInfo").innerText,
          selectedOverlay: p.data[p.data.length - 1],
          traceCount: p.data.length,
          tableSelected: Array.from(document.querySelectorAll("#regions tbody tr.selected")).map(el => el.innerText),
          tab: document.querySelector(".ram-main .nav-tabs li.active a").innerText,
          button: document.querySelector("#clearResidue").outerHTML
        };
      });
      console.error("Clear-selection diagnostics:", JSON.stringify(state));
      await page.screenshot({
        path: "benchmarks/output/ui-preview/selection-debug.png", fullPage: true
      });
      throw e;
    }
    // The NGL stage has a pick signal and a named highlight representation.
    await page.waitForFunction(() => window.getNGLStage &&
      window.getNGLStage("NGL") &&
      window.getNGLStage("NGL").getRepresentationsByName("ram-highlight")
        .list.length > 0, { timeout: 18000 });

    // Simulate an NGL hit on a known atom through its actual stage signal.
    // NGLVieweR's own handler expects getLabel(); our linked selector also
    // supports closestBondAtom for bond picking. No artificial Shiny input.
    const previousZoomCount = await page.evaluate(() => window.__ramZoomCalls.length);
    await page.evaluate(detail => {
      const stage = window.getNGLStage("NGL");
      stage.signals.clicked.dispatch({
        getLabel() { return "Selected backbone atom"; },
        closestBondAtom: {
          chainname: String(detail[0]), resno: Number(detail[1]),
          inscode: String(detail[2] || "")
        }
      });
    }, firstPoint);
    await page.waitForFunction(({ selected, minZoom }) => {
      const p = document.getElementById("plotly");
      return document.querySelector("#selectedResidueInfo strong") &&
        p.data[p.data.length - 1].x.length === 1 &&
        window.__ramStickSelections.includes(selected) &&
        window.__ramZoomCalls.length > minZoom &&
        window.__ramZoomCalls.some(call => call.sele === selected);
    }, { timeout: 20000 }, {
      selected: firstNglSelector, minZoom: previousZoomCount
    });

    // Check whether the browser delivered resize events and whether Plotly's
    // relayout handler actually received the new mobile dimensions.
    await page.evaluate(() => {
      window.__ramResizeEvents = 0;
      window.__ramRelayoutCalls = [];
      window.addEventListener("resize", () => window.__ramResizeEvents++);
      const original = window.Plotly.relayout;
      window.Plotly.relayout = function(node, update, ...rest) {
        if (node && node.id === "plotly") {
          window.__ramRelayoutCalls.push({ ...update });
        }
        return original.call(this, node, update, ...rest);
      };
    });
    await page.setViewport({ width: 390, height: 844, deviceScaleFactor: 1 });
    try {
      await page.waitForFunction(() => {
        const p = document.getElementById("plotly");
        if (!p || !p._fullLayout) return false;
        const width = p.clientWidth;
        const cardWidth = p.closest(".ram-chart-card").clientWidth;
        const targetHeight = Math.max(285, Math.min(690, Math.round(width + 35)));
        return width > 0 && width <= cardWidth + 2 &&
               Math.abs(p._fullLayout.width - width) <= 3 &&
               Math.abs(p._fullLayout.height - targetHeight) <= 3 &&
               p._fullLayout.xaxis._length >= 0.45 * width &&
               p._fullLayout.yaxis._length >= 0.45 * width;
      }, { timeout: 15000 });
    } catch (e) {
      const details = await page.evaluate(() => {
        const p = document.getElementById("plotly");
        return {
          clientWidth: p.clientWidth,
          styleHeight: p.style.height,
          computedHeight: getComputedStyle(p).height,
          plotClass: p.className,
          windowWidth: window.innerWidth,
          resizeEvents: window.__ramResizeEvents,
          relayoutCalls: window.__ramRelayoutCalls,
          activeTab: document.querySelector(".ram-main .nav-tabs li.active a").textContent,
          layoutWidth: p._fullLayout && p._fullLayout.width,
          layoutHeight: p._fullLayout && p._fullLayout.height,
          axisWidth: p._fullLayout && p._fullLayout.xaxis._length,
          axisHeight: p._fullLayout && p._fullLayout.yaxis._length,
          cardWidth: p.closest(".ram-chart-card").clientWidth
        };
      });
      console.error("Responsive layout diagnostics:", JSON.stringify(details));
      console.error("Browser page errors:", errors);
      await page.screenshot({
        path: "benchmarks/output/ui-preview/mobile-debug.png", fullPage: true
      });
      throw e;
    }
    const dimensions = await page.evaluate(() => ({
      viewport: document.documentElement.clientWidth,
      body: document.body.scrollWidth,
      main: document.querySelector(".ram-main").getBoundingClientRect().width,
      resultsTop: document.querySelector(".ram-main").getBoundingClientRect().top,
      settingsTop: document.querySelector(".ram-sidebar").getBoundingClientRect().top,
      settingsLinkVisible: document.querySelector(".ram-mobile-settings-link")
        .getBoundingClientRect().width > 0,
      viewerControlCount: document.querySelectorAll(".ram-viewer-options input[type=checkbox]").length
    }));
    assert.ok(dimensions.body <= dimensions.viewport + 3,
              "The mobile UI should not scroll horizontally.");
    assert.ok(dimensions.resultsTop < dimensions.settingsTop,
              "Analysis results must appear before settings on mobile.");
    assert.ok(dimensions.settingsLinkVisible,
              "The mobile results should link to analysis settings.");
    assert.equal(dimensions.viewerControlCount, 5,
                 "All 3D controls must stay beside the viewer.");
    await page.screenshot({
      path: "benchmarks/output/ui-preview/mobile-loaded.png", fullPage: true
    });
    if (errors.length) throw Error("Browser errors: " + errors.join(" | "));
    console.log(JSON.stringify({ desktop: plot, mobile: dimensions }));
  } finally {
    await browser.close();
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
