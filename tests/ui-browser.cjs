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
    // Use DOM collections directly; Shiny may replace individual controls
    // while a reactive update is being applied to the page.
    const cleanControls = await page.evaluate(() => {
      const styles = document.querySelectorAll('input[name="nglRepresentation"]');
      const layers = document.querySelectorAll(
        '.ram-viewer-options input[type="checkbox"]');
      return {
        compactInspector: document.querySelector(".ram-global-inspector")
          .classList.contains("is-empty"),
        previousHidden: getComputedStyle(document.getElementById("prevReview")).display
          === "none",
        reviewLabel: document.getElementById("nextReview").textContent.trim(),
        molecularStyles: Array.prototype.map.call(styles || [], item => item.value),
        selectedStyle: document.querySelector(
          'input[name="nglRepresentation"]:checked')?.value,
        visibleLayers: Array.prototype.every.call(layers || [], item =>
          item.closest(".ram-toggles")?.getBoundingClientRect().height > 0)
      };
    });
    assert.ok(cleanControls.compactInspector && cleanControls.previousHidden,
      "The empty residue inspector must not display unnecessary buttons.");
    assert.equal(cleanControls.reviewLabel,"Review issues");
    assert.deepEqual(cleanControls.molecularStyles,
      ["cartoon","ribbon","licorice","ball+stick","surface"]);
    assert.equal(cleanControls.selectedStyle, "cartoon",
      "Cartoon must be the default molecular representation.");
    assert.equal(cleanControls.visibleLayers, true,
      "All five molecular layer and motion switches must remain visible.");
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
    assert.equal(await page.$eval("#colorscheme", el => el.value),
                 "RamplotR", "RamplotR must be the publication-default palette.");
    assert.equal(await page.evaluate(() =>
      document.getElementById("plotly").data[0].fillcolor),
      "#D4ECE7", "The default contour must use the RamplotR identity palette.");
    // NGL loads and paints asynchronously after the Plotly response.
    // Give the viewer a moment to render before taking the desktop preview.
    await new Promise(resolve => setTimeout(resolve, 1600));
    await page.screenshot({
      path: "benchmarks/output/ui-preview/desktop-loaded.png", fullPage: true
    });
    // Display choices and secondary options are discoverable without
    // opening another disclosure or leaving the molecular view.
    assert.equal(await page.$eval('.ram-viewer-layers input[type="checkbox"]',
      nodes => nodes.length), 5,
      "The 3D viewer must expose all five layer and motion controls.");

    // A standard laptop should show the analysis, not a full-height landing
    // page. Focus mode is reversible and must preserve the live plot.
    await page.setViewport({ width: 1366, height: 768, deviceScaleFactor: 1 });
    await page.waitForFunction(() => {
      const app = document.querySelector(".ram-app");
      const p = document.querySelector("#plotly");
      return app && app.classList.contains("ram-has-data") &&
        p && p._fullLayout && Math.abs(p._fullLayout.width-p.clientWidth)<3;
    }, {timeout:15000});
    const laptop = await page.evaluate(() => {
      const intro = document.querySelector(".ram-intro");
      const source = document.querySelector(".ram-source");
      const p = document.querySelector("#plotly");
      return {
        introHidden: getComputedStyle(intro).display === "none",
        sourceHeight: source.getBoundingClientRect().height,
        plotBottom: p.getBoundingClientRect().bottom,
        viewerHeight: document.querySelector(".ram-ngl")
          .getBoundingClientRect().height,
        viewportHeight: window.innerHeight
      };
    });
    assert.ok(laptop.introHidden,
              "Loaded analysis should reclaim the introductory hero area.");
    assert.ok(laptop.sourceHeight < 100,
              "The loaded structure toolbar should remain compact.");
    assert.ok(laptop.plotBottom <= laptop.viewportHeight+95,
              "The laptop plot should fit mostly inside the first screen.");
    assert.ok(laptop.viewerHeight >= 290 && laptop.viewerHeight < 380,
              "The molecular viewer should fit a short laptop viewport.");
    await page.screenshot({
      path:"benchmarks/output/ui-preview/laptop-compact.png",fullPage:true
    });
    await page.click("#ram-toggle-settings");
    await page.waitForFunction(() => {
      const app = document.querySelector(".ram-app");
      const sidebar = document.querySelector(".ram-sidebar");
      return app.classList.contains("ram-focus-mode") &&
        getComputedStyle(sidebar).display === "none" &&
        document.querySelector("#ram-toggle-settings")
          .getAttribute("aria-expanded") === "false";
    });
    await page.screenshot({
      path:"benchmarks/output/ui-preview/laptop-focus.png",fullPage:true
    });
    await page.click("#ram-toggle-settings");
    await page.waitForFunction(() =>
      !document.querySelector(".ram-app").classList.contains("ram-focus-mode"));
    await page.setViewport({width:1440,height:940,deviceScaleFactor:1});
    await page.waitForFunction(() => {
      const p = document.querySelector("#plotly");
      return p._fullLayout && Math.abs(p._fullLayout.width-p.clientWidth)<3;
    }, {timeout:15000});

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
    await page.waitForFunction(() =>
      !document.querySelector(".ram-global-inspector").classList.contains("is-empty")
      && getComputedStyle(document.getElementById("prevReview")).display !== "none",
      {timeout:10000});
    await page.screenshot({
      path: "benchmarks/output/ui-preview/desktop-residue-zoom.png",
      fullPage: true
    });

    // The residue table should also select the same entry, and clicking
    // another row should move the plot highlight back to that residue.
    // Change the protein representation while retaining the orange residue
    // overlay; then restore cartoon before continuing linked-view tests.
    await page.evaluate(() =>
      document.querySelector('input[name="nglRepresentation"][value="licorice"]').click());
    await page.waitForFunction(() => {
      const group = window.getNGLStage("NGL")?.getRepresentationsByName("ram-chain-A");
      return group && group.list && group.list.some(item =>
        (item.repr?.type || "").toLowerCase() === "licorice");
    }, {timeout:20000});
    await page.evaluate(() =>
      document.querySelector('input[name="nglRepresentation"][value="cartoon"]').click());
    await page.waitForFunction(() => {
      const group = window.getNGLStage("NGL")?.getRepresentationsByName("ram-chain-A");
      const highlight = window.getNGLStage("NGL")?.getRepresentationsByName("ram-highlight");
      return group && group.list && group.list.some(item =>
        (item.repr?.type || "").toLowerCase() === "cartoon") &&
        highlight && highlight.list.length > 0;
    }, {timeout:20000});

    await page.click('.nav-tabs a[data-value="residues"]');
    await page.waitForSelector("#regions table tbody tr", { timeout: 18000 });
    await page.waitForFunction(() =>
      document.querySelector("#regions tbody tr.selected"), { timeout: 18000 });
    // One physical table (no DataTables scroll-head clone) must align the
    // column headings with the corresponding residue values.
    const geometry = await page.evaluate(() => {
      const table = document.querySelector("#regions table");
      const head = table.querySelectorAll("thead th");
      const cells = table.querySelectorAll("tbody tr:first-child td");
      return [0,1,2,3,4,5,6,7].map(i => ({
        heading: head[i].getBoundingClientRect().left,
        value: cells[i].getBoundingClientRect().left
      }));
    });
    assert.ok(geometry.every(col => Math.abs(col.heading - col.value) <= 4),
              "Residue table header and body columns must align.");
    await page.screenshot({
      path:"benchmarks/output/ui-preview/residue-table.png",fullPage:true
    });
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

    // Sequence navigation belongs to the plot tab. The short all-chain
    // overview stays visible while the letter-level navigator is collapsed.
    await page.waitForSelector("#ram-sequence-panel > summary");
    const shortView = await page.evaluate(() => ({
      overviewChains: document.querySelectorAll(".ram-sequence-overview-chain").length,
      collapsed: !document.getElementById("ram-sequence-panel").open,
      plotTab: document.querySelector('.nav-tabs li.active a').getAttribute("data-value")
    }));
    assert.equal(shortView.overviewChains, 1,
                 "The small single-chain 1CRN fixture has one overview row.");
    assert.ok(shortView.collapsed && shortView.plotTab === "plot",
              "Sequence must be integrated into the plot tab and initially collapsed.");
    await page.click("#ram-sequence-panel > summary");
    await page.waitForSelector("#sequenceView .ram-seq-res", {timeout:18000});
    const sequencePick = await page.evaluate(first => {
      const buttons = Array.from(document.querySelectorAll(".ram-seq-res"));
      const other = buttons.find(b =>
        b.dataset.chain !== String(first[0]) ||
        Number(b.dataset.resi) !== Number(first[1]) ||
        b.dataset.insertion !== String(first[2] || ""));
      if (!other) return null;
      other.click();
      return other.dataset.resi;
    }, firstPoint);
    assert.ok(sequencePick, "Sequence navigator needs multiple residues.");
    await page.waitForFunction(resi => {
      const chosen = document.querySelector("#selectedResidueInfo strong");
      return chosen && chosen.textContent.includes(" " + resi);
    }, {timeout:15000}, sequencePick);
    assert.equal(await page.$eval('.nav-tabs li.active a',
      el => el.getAttribute("data-value")), "plot",
      "Picking a residue must preserve the visible plot and NGL viewer.");
    await page.screenshot({
      path:"benchmarks/output/ui-preview/sequence-integrated.png", fullPage:true
    });
    await page.click("#ram-sequence-panel > summary");

    // Comparing 1CRN with itself tests sequence alignment without relying
    // on any additional external downloads inside the Shiny session.
    await page.click('.nav-tabs a[data-value="compare"]');
    await page.evaluate(() =>
      document.querySelector('input[name="compareInputSource"][value="upload"]').click());
    await page.waitForFunction(() => {
      const element = document.getElementById("ram-compare-upload");
      return element && !element.classList.contains("is-hidden");
    }, {timeout:15000});
    await (await page.$("#compareFile")).uploadFile(
      path.resolve("benchmarks/output/ui-preview/1CRN.pdb"));
    await new Promise(resolve => setTimeout(resolve, 1500));
    await page.click("#compareSubmit");
    await page.waitForSelector("#comparison tbody tr", {timeout:25000});
    await page.waitForFunction(() => {
      const p = document.getElementById("comparePlot");
      return p && p.data && p.data.length >= 2;
    }, {timeout:18000});
    await page.screenshot({
      path:"benchmarks/output/ui-preview/compare-self.png",fullPage:true
    });
    await page.click('.nav-tabs a[data-value="plot"]');

    // Multi-chain experimental fixture: the compact overview and the
    // expanded residue navigator must each show A, B, C and D together.
    await (await page.$("#structfile")).uploadFile(
      path.resolve("benchmarks/output/ui-preview/1BBB.pdb"));
    await new Promise(resolve => setTimeout(resolve, 1500));
    await page.click("#submit");
    await page.waitForFunction(() => {
      const badge = document.getElementById("ram-current-structure");
      const rows = document.querySelectorAll(".ram-sequence-overview-chain");
      return badge && badge.textContent.includes("1BBB") && rows.length === 4;
    }, {timeout:45000});
    // Chain selection is updated asynchronously after replacing a
    // single-chain structure. Do not mistake a transient Chain-A-only plot
    // for a completed four-chain analysis.
    await page.waitForFunction(() => {
      const plot = document.getElementById("plotly");
      const badge = document.getElementById("ram-current-structure");
      const chainTraces = (plot?.data || [])
        .filter(trace => /^Chain [A-D]$/.test(trace.name || ""));
      const stage = typeof window.getNGLStage === "function"
        ? window.getNGLStage("NGL") : null;
      const chainD = stage && stage.getRepresentationsByName("ram-chain-D");
      return badge && badge.textContent.includes("1BBB") &&
        chainTraces.length === 4 &&
        chainD && chainD.list && chainD.list.length > 0;
    }, {timeout:30000});
    // Let NGL's WebGL renderer complete a frame before archiving screenshots.
    await new Promise(resolve => setTimeout(resolve, 900));
    const overviewNames = await page.$eval(
      ".ram-sequence-overview-chain .ram-sequence-chain-name",
      nodes => nodes.map(n => n.textContent.trim()));
    assert.deepEqual(overviewNames, ["Chain A", "Chain B", "Chain C", "Chain D"],
                     "The collapsed navigator should show every protein chain.");
    await page.click("#ram-sequence-panel > summary");
    await page.waitForFunction(() =>
      document.querySelectorAll("#sequenceView .ram-sequence-chain").length === 4,
      {timeout:15000});
    const allChains = await page.evaluate(() =>
      Array.from(document.querySelectorAll(
        "#sequenceView .ram-sequence-chain")).map(group => ({
          title: group.querySelector("strong").textContent,
          residues: group.querySelectorAll(".ram-seq-res").length
        })));
    assert.equal(allChains.length, 4);
    assert.ok(allChains.every(group => group.residues > 0),
              "Every chain should provide clickable residue navigation.");
    await page.screenshot({
      path:"benchmarks/output/ui-preview/all-chains-expanded.png", fullPage:true
    });
    await page.click("#sequenceView .ram-sequence-chain:last-child .ram-seq-res:nth-child(2)");
    await page.waitForFunction(() => {
      const selected = document.querySelector("#selectedResidueInfo strong");
      return selected && selected.textContent.includes("Chain D");
    }, {timeout:15000});
    await page.click("#ram-sequence-panel > summary");

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
