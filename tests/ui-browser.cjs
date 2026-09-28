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
    page.on("pageerror", error => errors.push(error.message));
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

    await page.setViewport({ width: 390, height: 844, deviceScaleFactor: 1 });
    await page.waitForFunction(() => {
      const p = document.getElementById("plotly");
      if (!p || !p._fullLayout) return false;
      const width = p.clientWidth;
      const targetHeight = Math.max(285, Math.min(690, Math.round(width + 35)));
      return Math.abs(p._fullLayout.height - targetHeight) <= 3 &&
             p._fullLayout.xaxis._length >= 0.45 * width &&
             p._fullLayout.yaxis._length >= 0.45 * width;
    }, { timeout: 15000 });
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
