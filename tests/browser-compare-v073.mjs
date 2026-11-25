// Fixed-camera original v0.72 vs common CPU reference. Requires local server.
import { chromium } from "playwright";
import { readFileSync, writeFileSync, mkdirSync, existsSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { join } from "node:path";
const out = "test-results/comparison-v073",
  cache = process.env.ATLAS_CACHE || "/tmp/upper-limb-atlas";
mkdirSync(out, { recursive: true });
const browser = await chromium.launch({
  executablePath:
    process.env.CHROME_PATH ||
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
  headless: true,
  args: [
    "--ignore-certificate-errors",
    "--use-angle=swiftshader",
    "--enable-unsafe-swiftshader",
  ],
});
const reference = execFileSync(
    "git",
    ["show", "archive/v0.72.0-before-common-engine:index.html"],
    { encoding: "utf8" },
  ),
  current = readFileSync("index.html", "utf8");
const cameras = new Map(),
  results = [];
try {
  for (const [version, source] of [
    ["v072", reference],
    ["v073", current],
  ]) {
    const page = await browser.newPage({
        viewport: { width: 1280, height: 900 },
      }),
      errors = [];
    page.on("pageerror", (e) => errors.push(String(e)));
    await page.route("**/*.stl", async (route) => {
      const file = new URL(route.request().url()).pathname.split("/").at(-1),
        p = join(cache, file);
      if (!existsSync(p))
        throw Error("Populate atlas cache with browser-v073 first: " + file);
      await route.fulfill({
        body: readFileSync(p),
        contentType: "application/octet-stream",
        headers: { "access-control-allow-origin": "*" },
      });
    });
    const instrument = String.raw`window.__compare={set({camera,elev,plane,yaw}){$('elev').value=elev;$('plane').value=plane;for(const spec of V60_ALL)$(spec.id).checked=spec.id.startsWith('v60Pec')&&['v60PecClav','v60PecStern','v60PecAbd'].includes(spec.id)||spec.deltoid;for(const id of ['muscles','scapularMuscles','brachioradialis','pronatorTeres','supinator','flexorCarpiRadialis','flexorCarpiUlnaris','palmarisLongus','flexorDigitorumSuperficialis','flexorDigitorumProfundus',...V43_MUSCLES])$(id).checked=false;update();if($('v73HeadVisible'))$('v73HeadVisible').checked=false;if($('v73NeckMusclesVisible'))$('v73NeckMusclesVisible').checked=false;if(camera)Object.assign(state,camera);else {state.yaw=yaw;state.pitch=.1;fitPose()}request();return {yaw:state.yaw,pitch:state.pitch,zoom:state.zoom,target:state.target}},draw(){render()}};`;
    const html = source
      .replace('<option value="cpu">', '<option value="cpu" selected>')
      .replace(/\}\)\(\);\s*<\/script>/, () => instrument + "})();</script>");
    await page.route("**/index.html", (r) =>
      r.fulfill({ body: html, contentType: "text/html" }),
    );
    await page.goto(process.env.QA_URL || "http://localhost:8730/index.html");
    await page.waitForFunction(
      () =>
        window.__labDeepMuscles?.().loaded.length === 17 &&
        (!window.__labCervicalV073 ||
          window.__labCervicalV073().neckMuscles === 6),
      null,
      { timeout: 240000 },
    );
    await page
      .locator("#v63LoadShow")
      .click({ timeout: 1000 })
      .catch(() => {});
    for (const [elev, plane] of [
      [0, 0],
      [150, 90],
    ])
      for (const [view, yaw] of [
        ["front", Math.PI],
        ["side", Math.PI / 2],
        ["back", 0],
      ]) {
        const key = [elev, plane, view].join("-"),
          poseCamera = await page.evaluate((c) => window.__compare.set(c), {
            camera: cameras.get(key),
            elev,
            plane,
            yaw,
          });
        if (!cameras.has(key)) cameras.set(key, poseCamera);
        if (JSON.stringify(poseCamera) !== JSON.stringify(cameras.get(key)))
          throw Error("Camera mismatch");
        await page.evaluate(() => window.__compare.draw());
        const rect = await page.locator("#viewer").boundingBox();
        await page.screenshot({
          clip: rect,
          path: join(out, `${version}-${key}.png`),
          timeout: 120000,
        });
        results.push({
          version,
          focus:
            "pectoralis and deltoid subdivisions; other muscles hidden equally",
          view,
          camera: poseCamera,
          pose: { right: { elev, plane }, left: { elev: 0, plane: 0 } },
          renderer:
            version === "v073"
              ? "CPU reference, WebGL rasterization"
              : "original archived v0.72 GPU fallback",
          errors,
        });
        writeFileSync(
          join(out, "comparison.json"),
          JSON.stringify(results, null, 2),
        );
        console.log(version, key);
      }
    await page.close();
  }
  writeFileSync(join(out, "comparison.json"), JSON.stringify(results, null, 2));
  console.log(results);
} finally {
  await browser.close();
}
