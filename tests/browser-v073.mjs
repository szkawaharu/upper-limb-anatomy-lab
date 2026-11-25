// Actual-atlas desktop Chromium QA. Optional network/cache, no phone claims.
// npm ci; npm run serve; CHROME_PATH=... npm run test:browser
import { chromium } from "playwright";
import { readFileSync, writeFileSync, mkdirSync, existsSync } from "node:fs";
import { join } from "node:path";
import assert from "node:assert/strict";
import { trackStep } from "./browser-step-diagnostics.mjs";
const numericOnly = process.env.QA_NUMERIC_ONLY === "1";
const out = process.env.QA_OUTPUT || "test-results/browser-v073",
  cache = process.env.ATLAS_CACHE || "/tmp/upper-limb-atlas";
mkdirSync(out, { recursive: true });
mkdirSync(cache, { recursive: true });
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
const page = await browser.newPage({ viewport: { width: 1280, height: 900 } }),
  errors = [],
  warnings = [];
page.on("pageerror", (e) => errors.push(String(e)));
page.on("console", (m) => {
  if (m.type() === "warning") warnings.push(m.text());
});
await page.route("**/*.stl", async (route) => {
  const file = new URL(route.request().url()).pathname.split("/").at(-1),
    p = join(cache, file);
  try {
    if (!existsSync(p)) {
      const r = await fetch(
        "https://raw.githubusercontent.com/Kevin-Mattheus-Moerman/BodyParts3D/2c83c698788cff57f6317bddcf712e7acfbea7d8/assets/BodyParts3D_data/stl/" +
          file,
      );
      if (!r.ok) throw Error(r.status);
      writeFileSync(p, Buffer.from(await r.arrayBuffer()));
    }
    await route.fulfill({
      body: readFileSync(p),
      contentType: "application/octet-stream",
      headers: { "access-control-allow-origin": "*" },
    });
  } catch (e) {
    await route.abort();
  }
});
const instrument = String.raw`
let qaSkipDraw=false;const qaOriginalRender=render;render=()=>{if(!qaSkipDraw)qaOriginalRender();else{lastW=canvas.clientWidth;lastH=canvas.clientHeight;proj=ortho(state.zoom,lastW/lastH,.01,100);const eye=add(state.target,[18*Math.sin(state.yaw)*Math.cos(state.pitch),18*Math.sin(state.pitch),18*Math.cos(state.yaw)*Math.cos(state.pitch)]);view=lookAt(eye,state.target);pv=mm(proj,view)}};
const qaExport=(f,id)=>{const B=f.binding,samples=[];for(let i=0;i<B.weights.length;i+=Math.max(1,Math.floor(B.weights.length/80))){const p=Array.from(B.rest.slice(i*3,i*3+3)),t=B.weights[i],l=B.lane[i],materialT=B.material[i],q=MuscleDeformation.point(f,p,t,l,materialT);samples.push({vertexId:i,p,t,l,materialT,q,originBlend:B.originBlend?.[i]||0,originAnchor:f.originAt?f.originAt(p):null})}return {id,name:B.definition.name||id,side:B.definition.side||"right",samples,state:{a:f.a,b:f.b,rest:Array.from(f.restFlat),now:Array.from(f.nowFlat),restControls:B.field.controls,posedControls:f.posed,sheet:B.sheet,blendRange:B.blendRange,lanes:B.field.controls.length,tendonStart:B.definition.tissues.tendonStart,reference:B.reference,contact:f.contact,multiBone:!!f.originAt,originMix:f.originMix}}};
window.__qa={
 ensureGPU(){if(!v73Renderer)v73Renderer=new MuscleRenderer(gl,shaderFS,(v,f)=>gpu(mesh(v,f)))},
 skipRender(value){qaSkipDraw=value;if(!value)request()},
 pose(p){for(const [k,v] of Object.entries(p))$(k).value=v;update()},
 both(right,left){v56ChooseSide('right');this.pose(right);v56ChooseSide('left');this.pose(left);v56ChooseSide('right')},
 neck(n){v73SetNeck(n)},
 view(yaw){state.yaw=yaw;state.pitch=.1;fitPose();request()},
 boneTrace(elev,plane,axial=0){const data=v61GuideData(),raw=shoulderPoseFor(elev,plane,axial,true,true),final=v72ShoulderGuide(raw.qC,raw.qS,raw.AC,elev,plane,raw.qH,axial,data);const compact=g=>({qC:g?.qC,qS:g?.qS,AC:g?.AC,params:g?.params,blend:g?.v67?.blend,reason:g?.v72?.reason,contactParams:g?.v72?.params,before:g?.v72?.beforeScore,after:g?.v72?.afterScore,regions:g?.v72?.after?.regions||g?.metrics?.regions});return {elev,plane,axial,raw:compact(raw),v66:compact(v66Last),v67:compact(v67Last),v72:compact(final)}},
 camera(yaw, prior){if(prior)Object.assign(state,prior);else{state.yaw=yaw;state.pitch=.1;fitPose()}return {yaw:state.yaw,pitch:state.pitch,zoom:state.zoom,target:state.target.slice()}},
 sample(){return V60_ALL.map(spec=>{const m=v58Meshes.get(spec.id);return m?.common?qaExport(v73MuscleState(spec,m),spec.id):{id:spec.id,error:'missing'}})},
 neckSamples(){return [...v73NeckMuscles.values()].map(m=>qaExport(v73NeckMuscleState(m),m.spec.id))},

 audit(ids){return V60_ALL.filter(s=>ids.includes(s.id)).map(s=>{const m=v58Meshes.get(s.id),f=v73MuscleState(s,m);return DeformationAudit.audit(f,MuscleDeformation.deform(f))})},
 headPoint(){const pts=v73HeadBounds().map(projectForPicking).filter(Boolean);if(!pts.length)return null;const r=canvas.getBoundingClientRect();return [r.left+(Math.min(...pts.map(p=>p[0]))+Math.max(...pts.map(p=>p[0])))/2,r.top+(Math.min(...pts.map(p=>p[1]))+Math.max(...pts.map(p=>p[1])))/2]},
 gpuSource(){return MuscleGPU.vertexShader()},
 status(){return {loaded:window.__labDeepMuscles(),head:window.__labCervicalV073(),arms:window.__labUpperBody(),camera:{yaw:state.yaw,pitch:state.pitch},backend:v73Renderer?.backend,reason:v73Renderer?.reason,glError:gl.getError(),render:window.__labRefinementV0601().renderCPU}},
 render(){render()},
};
`;
await page.route("**/index.html", async (route) => {
  const html = readFileSync("index.html", "utf8").replace(
    /\}\)\(\);\s*<\/script>/,
    () =>
      instrument.replace(
        "let qaSkipDraw=false",
        "let qaSkipDraw=" + numericOnly,
      ) + "})();</script>",
  );
  await route.fulfill({ body: html, contentType: "text/html" });
});
try {
  await page.goto(process.env.QA_URL || "http://localhost:8730/index.html");
  await page.waitForFunction(
    () =>
      window.__labDeepMuscles?.().loaded.length === 17 &&
      window.__labCervicalV073?.().skullParts === 22 &&
      window.__labCervicalV073?.().neckMuscles === 6,
    null,
    { timeout: 240000 },
  );
  await page
    .locator("#v63LoadShow")
    .click({ timeout: 1000 })
    .catch(() => {});
  console.log("Actual atlas loaded");
  await page.evaluate(() => window.__qa.ensureGPU());
  const initial = await page.evaluate(() => window.__qa.status());
  assert.equal(initial.backend, "gpu");
  assert.equal(initial.glError, 0);
  assert.deepEqual(
    (await page.evaluate(() => window.__qa.sample())).filter((x) => x.error),
    [],
    "missing common binding",
  );
  if(process.env.QA_TRACE_ONLY==="1") {
    const trace=[];
    for(let elev=114;elev<=129;elev++)trace.push(await page.evaluate(e=>window.__qa.boneTrace(e,0),elev));
    writeFileSync(join(out,"bone-trace.json"),JSON.stringify(trace,null,2));
    for(const t of trace)console.log("BONE_TRACE "+JSON.stringify(t));
    const cameras={};
    await page.evaluate(()=>window.__qa.skipRender(false));
    for(const elev of [120,123]){
      await page.evaluate(e=>window.__qa.both({elev:e,plane:0},{elev:0,plane:0}),elev);
      for(const [name,yaw] of [["back",0],["side",Math.PI/2],["front",Math.PI]]){
        const camera=await page.evaluate(([y,c])=>window.__qa.camera(y,c),[yaw,cameras[name]||null]);
        if(cameras[name])assert.deepEqual(camera,cameras[name]);else cameras[name]=camera;
        await page.evaluate(()=>window.__qa.render());
        const rect=await page.locator("#viewer").boundingBox();
        await page.screenshot({clip:rect,path:join(out,`baseline-${elev}-0-${name}.png`),timeout:120000});
      }
    }
    writeFileSync(join(out,"bone-trace-camera.json"),JSON.stringify(cameras,null,2));
    assert.equal(errors.length,0,errors.join("\n"));
  } else {
  const rows = [],
    paths = [],
    cases = [
      [0, 0],
      [30, 0],
      [60, 0],
      [90, 0],
      [120, 0],
      [150, 0],
      [180, 0],
      [30, 90],
      [60, 90],
      [90, 90],
      [120, 90],
      [150, 90],
      [180, 90],
    ];
  for (const [elev, plane] of cases) {
    console.log("Pose", elev, plane);
    await page.evaluate(
      ([e, p]) =>
        window.__qa.both({ elev: e, plane: p }, { elev: e, plane: p }),
      [elev, plane],
    );
    const samples = await page.evaluate(() => window.__qa.sample());
    paths.push({ elev, plane, samples });
    if ([0, 90, 150, 180].includes(elev)) {
      rows.push({
        elev,
        plane,
        audit: await page.evaluate(() =>
          window.__qa.audit([
            "v60PecClav",
            "v60PecStern",
            "v60PecAbd",
            "v60DeltAnt",
            "v60DeltMid",
            "v60DeltPost",
          ]),
        ),
      });
      if (!numericOnly)
        for (const [name, yaw] of [
          ["front", Math.PI],
          ["side", Math.PI / 2],
          ["back", 0],
          ["oblique", Math.PI * 0.78],
        ]) {
          await page.evaluate((y) => window.__qa.view(y), yaw);
          await page.waitForTimeout(80);
          await page.screenshot({
            path: join(out, `${elev}-${plane}-${name}.png`),
          });
        }
    }
  }
  // Continuous forward/abduction sweeps and exact inverse return; samples are
  // material vertices, not refitted geometry or a pose-history-dependent solver.
  writeFileSync(join(out, "pose-audit.json"), JSON.stringify(rows, null, 2));
  await page.evaluate(() => window.__qa.skipRender(true));
  let maxStepMm = 0,
    maxReturnMm = 0;
  const maxima = { forward: {}, reverse: {} };
  for (const plane of [0, 90]) {
    console.log("Sweep", plane);
    let prior = null;
    const forward = [];
    for (let elev = 0; elev <= 180; elev += 3) {
      await page.evaluate((p) => window.__qa.pose(p), { elev, plane });
      const samples = await page.evaluate(() => window.__qa.sample());
      for (let m = 0; m < samples.length; m++)
        for (let i = 0; i < samples[m].samples.length; i++) {
          const q = samples[m].samples[i].q;
          assert(q.every(Number.isFinite));
          if (prior) {
            const d = trackStep(maxima.forward, prior[m], samples[m], i, {plane, fromElev:elev-3, toElev:elev, direction:"forward"});
            maxStepMm = Math.max(maxStepMm, d);
          }
        }
      forward.push(samples);
      prior = samples;
    }
    prior = null;
    for (let elev = 180; elev >= 0; elev -= 3) {
      await page.evaluate((p) => window.__qa.pose(p), { elev, plane });
      const samples = await page.evaluate(() => window.__qa.sample()),
        ref = forward[elev / 3];
      for (let m = 0; m < samples.length; m++)
        for (let i = 0; i < samples[m].samples.length; i++) {
          if (prior) trackStep(maxima.reverse, prior[m], samples[m], i, {plane, fromElev:elev+3, toElev:elev, direction:"reverse"});
          maxReturnMm = Math.max(
            maxReturnMm,
            Math.hypot(
              ...samples[m].samples[i].q.map(
                (x, k) => x - ref[m].samples[i].q[k],
              ),
            ) / 0.00555,
          );
        }
      prior = samples;
    }
  }
  const stepDiagnostics = { samplePolicy: "original stride floor(vertexCount/80), minimum 1; zero-based clustered mesh vertex IDs", scale: 0.00555, forward: Object.values(maxima.forward).sort((a,b)=>b.mm-a.mm), reverse: Object.values(maxima.reverse).sort((a,b)=>b.mm-a.mm) };
  writeFileSync(join(out,"step-diagnostics.json"),JSON.stringify(stepDiagnostics,null,2));
  console.log("Largest forward and reverse sampled steps", JSON.stringify({forward:stepDiagnostics.forward[0],reverse:stepDiagnostics.reverse[0]}));
  await page.evaluate((v) => window.__qa.skipRender(v), numericOnly);
  console.log("Neck gestures");
  const neck = [];
  await page.evaluate(() =>
    window.__qa.both(
      { elev: 120, plane: 90, elbow: 70 },
      { elev: 90, plane: 20, elbow: 100 },
    ),
  );
  const arms = await page.evaluate(() => window.__labUpperBody().poses);
  await page.evaluate(() => window.__qa.skipRender(true));
  for (const angle of [-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75]) {
    await page.evaluate((a) => window.__qa.neck(a), angle);
    assert.deepEqual(
      await page.evaluate(() => window.__labUpperBody().poses),
      arms,
    );
    neck.push(await page.evaluate(() => window.__labCervicalV073()));
  }
  await page.evaluate((v) => window.__qa.skipRender(v), numericOnly);
  const drag = [];
  for (const yaw of [0, Math.PI / 2, Math.PI]) {
    await page.evaluate((y) => {
      window.__qa.neck(0);
      window.__qa.view(y);
    }, yaw);
    await page.waitForTimeout(50);
    const xy = await page.evaluate(() => window.__qa.headPoint());
    assert(xy);
    const before = await page.evaluate(() => window.__qa.status());
    await page.mouse.move(...xy);
    await page.mouse.down();
    await page.mouse.move(xy[0] + 75, xy[1], { steps: 8 });
    await page.mouse.up();
    const after = await page.evaluate(() => window.__qa.status());
    assert(Math.abs(after.head.angle) > 0);
    assert.deepEqual(after.arms.poses, before.arms.poses);
    assert.deepEqual(after.camera, before.camera);
    drag.push({ yaw, angle: after.head.angle });
  }
  // WebGL2 transform feedback executes the production GLSL mapping and reads
  // float positions back; independent CPU results above are the reference.
  console.log("GPU readback");
  paths.push({
    neck: await page.evaluate(() => window.__labCervicalV073().angle),
    samples: await page.evaluate(() => window.__qa.neckSamples()),
  });
  const gpu = await page.evaluate(
    ({ source, sets }) => {
      const canvas = document.createElement("canvas"),
        g = canvas.getContext("webgl2");
      if (!g) return { available: false };
      let maxMm = 0,
        count = 0;
      const compile = (type, s) => {
        const x = g.createShader(type);
        g.shaderSource(x, s);
        g.compileShader(x);
        if (!g.getShaderParameter(x, g.COMPILE_STATUS))
          throw Error(g.getShaderInfoLog(x));
        return x;
      };
      const vs =
          "#version 300 es\n" +
          source.replaceAll("attribute ", "in ").replaceAll("varying ", "out "),
        fs =
          "#version 300 es\nprecision highp float;out vec4 color;void main(){color=vec4(1.);}";
      const p = g.createProgram();
      g.attachShader(p, compile(g.VERTEX_SHADER, vs));
      g.attachShader(p, compile(g.FRAGMENT_SHADER, fs));
      g.transformFeedbackVaryings(p, ["W"], g.INTERLEAVED_ATTRIBS);
      g.linkProgram(p);
      if (!g.getProgramParameter(p, g.LINK_STATUS))
        throw Error(g.getProgramInfoLog(p));
      g.useProgram(p);
      const u = (n) => g.getUniformLocation(p, n),
        I = new Float32Array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]);
      const mat = (a) => {
        const [x, y, z, w] = a.q,
          rot = [
            1 - 2 * y * y - 2 * z * z,
            2 * x * y + 2 * z * w,
            2 * x * z - 2 * y * w,
            0,
            2 * x * y - 2 * z * w,
            1 - 2 * x * x - 2 * z * z,
            2 * y * z + 2 * x * w,
            0,
            2 * x * z + 2 * y * w,
            2 * y * z - 2 * x * w,
            1 - 2 * x * x - 2 * y * y,
            0,
            0,
            0,
            0,
            1,
          ];
        for (let k = 0; k < 3; k++)
          rot[12 + k] =
            a.p[k] -
            rot[k] * a.rest[0] -
            rot[4 + k] * a.rest[1] -
            rot[8 + k] * a.rest[2];
        return rot;
      };
      g.uniformMatrix4fv(u("model"), false, I);
      g.uniformMatrix4fv(u("mvp"), false, I);
      for (const set of sets)
        for (const row of set.samples) {
          if (row.state.multiBone && !row.state.originMix) continue;
          const s = row.state,
            n = row.samples.length,
            buffers = [];
          for (const [name, data, size] of [
            ["pos", row.samples.flatMap((x) => x.p), 3],
            ["normal", row.samples.flatMap(() => [0, 0, 1]), 3],
            ["weight", row.samples.map((x) => x.t), 1],
            ["lane", row.samples.map((x) => x.l), 1],
            ["materialT", row.samples.map((x) => x.materialT), 1],
            ["originBlend", row.samples.map((x) => x.originBlend || 0), 1],
          ]) {
            const a = g.getAttribLocation(p, name);
            if (a < 0) continue;
            const b = g.createBuffer();
            buffers.push(b);
            g.bindBuffer(g.ARRAY_BUFFER, b);
            g.bufferData(g.ARRAY_BUFFER, new Float32Array(data), g.STATIC_DRAW);
            g.enableVertexAttribArray(a);
            g.vertexAttribPointer(a, size, g.FLOAT, false, 0, 0);
          }
          for (const [name, data] of [
            ["restControls[0]", s.rest],
            ["nowControls[0]", s.now],
          ]) {
            const a = new Float32Array(108);
            a.set(data);
            g.uniform3fv(u(name), a);
          }
          g.uniformMatrix4fv(u("anchorA"), false, mat(s.a));
          g.uniformMatrix4fv(u("anchorB"), false, mat(s.b));
          g.uniform2fv(u("blendRange"), s.blendRange);
          g.uniform1f(u("originMixOn"), s.originMix ? 1 : 0);
          g.uniform4fv(u("originQA"), s.originMix?.q || s.a.q);
          g.uniformMatrix4fv(u("originAnchor"), false, mat(s.originMix || s.a));
          g.uniform4fv(u("qa"), s.a.q);
          g.uniform4fv(u("qb"), s.b.q);
          g.uniform1f(u("lanes"), s.lanes);
          g.uniform1f(u("sheet"), s.sheet ? 1 : 0);
          g.uniform1f(u("tendonStart"), s.tendonStart);
          g.uniform3fv(u("reference"), s.reference);
          const c = s.contact;
          g.uniform1f(u("contactOn"), c ? 1 : 0);
          if (c) {
            for (const [name, v] of Object.entries({
              contactProximal: c.proximal ? 1 : 0,
              unitScale: c.scale,
              contactIn: c.inMm ?? 2,
              contactOut: c.outMm ?? 6,
              contactMax: c.maxMm ?? 20,
            }))
              g.uniform1f(u(name), v);
            g.uniform3fv(u("contactCenter"), c.center);
            g.uniform3fv(u("contactRestCenter"), c.restCenter || c.center);
            g.uniform3fv(u("contactRadii"), c.radii);
          }
          const b = g.createBuffer();
          g.bindBuffer(g.TRANSFORM_FEEDBACK_BUFFER, b);
          g.bufferData(g.TRANSFORM_FEEDBACK_BUFFER, n * 12, g.STREAM_READ);
          g.bindBufferBase(g.TRANSFORM_FEEDBACK_BUFFER, 0, b);
          g.enable(g.RASTERIZER_DISCARD);
          g.beginTransformFeedback(g.POINTS);
          g.drawArrays(g.POINTS, 0, n);
          g.endTransformFeedback();
          g.disable(g.RASTERIZER_DISCARD);
          const actual = new Float32Array(n * 3);
          g.getBufferSubData(g.TRANSFORM_FEEDBACK_BUFFER, 0, actual);
          for (let i = 0; i < n; i++) {
            const q = row.samples[i].q;
            maxMm = Math.max(
              maxMm,
              Math.hypot(...q.map((x, k) => x - actual[i * 3 + k])) / 0.00555,
            );
            count++;
          }
          g.bindBufferBase(g.TRANSFORM_FEEDBACK_BUFFER, 0, null);
          g.deleteBuffer(b);
          buffers.forEach((b) => g.deleteBuffer(b));
        }
      return { available: true, count, maxMm, error: g.getError() };
    },
    { source: await page.evaluate(() => window.__qa.gpuSource()), sets: paths },
  );
  assert(
    gpu.available && gpu.maxMm < 1 && gpu.error === 0,
    "CPU/GPU mismatch " + JSON.stringify(gpu),
  );
  assert(maxReturnMm < 0.001, "return drift");
  assert.equal(errors.length, 0, errors.join("\n"));
  writeFileSync(
    join(out, "report.json"),
    JSON.stringify(
      {
        initial,
        rows,
        neck,
        drag,
        gpu,
        maxStepMm,
        maxReturnMm,
        stepDiagnostics,
        errors,
        warnings,
        numericOnly,
        device: "desktop Chromium / SwiftShader; physical smartphone untested",
      },
      null,
      2,
    ),
  );
  console.log(JSON.stringify({ gpu, maxStepMm, maxReturnMm, errors }));
  }
} catch (error) {
  writeFileSync(
    join(out, "failure.json"),
    JSON.stringify(
      { error: String(error), stack: error.stack, errors, warnings },
      null,
      2,
    ),
  );
  throw error;
} finally {
  await browser.close();
}
