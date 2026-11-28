// Actual-atlas desktop Chromium QA. Optional network/cache, no phone claims.
// npm ci; npm run serve; CHROME_PATH=... npm run test:browser
import { chromium } from "playwright";
import { readFileSync, writeFileSync, mkdirSync, existsSync, readdirSync } from "node:fs";
import { createHash } from "node:crypto";
import { join } from "node:path";
import assert from "node:assert/strict";
import { staticServer } from "./helpers/static-server.mjs";
const numericOnly = process.env.QA_NUMERIC_ONLY === "1";
const levatorExperiment=process.env.QA_LEVATOR_EXPERIMENT==="1";
const levatorMaskFile=process.env.QA_LEVATOR_MASK_FILE;
const levatorMasks=levatorMaskFile?JSON.parse(readFileSync(levatorMaskFile,"utf8")):null;
const shoulderCompare=process.env.QA_SHOULDER_MUSCLE_COMPARE==="1"||levatorExperiment;
const shoulderPreview=process.env.QA_SHOULDER_PREVIEW==="1";
const previewFile=process.env.QA_PREVIEW_FILE||"test-results/v076-blend-preview/field.json";
const preview=shoulderPreview?JSON.parse(readFileSync(previewFile,"utf8")):null;
if(shoulderPreview){
 assert(shoulderCompare,"Preview is confined to the explicit comparison harness");
 assert(!process.env.QA_URL,"Preview must run against the local checkout");
 assert(preview.maxRigidErrorBoundMm<1&&preview.productAccepted===false);
 for(const [file,hash] of Object.entries(preview.sourceHashes))assert.equal(createHash('sha256').update(readFileSync(file)).digest('hex'),hash,'Preview input mismatch: '+file);
}
const out = process.env.QA_OUTPUT || "test-results/browser-v073",
  cache = process.env.ATLAS_CACHE || "/tmp/upper-limb-atlas";
mkdirSync(out, { recursive: true });
mkdirSync(cache, { recursive: true });
const sourceFiles=['index.html',...readdirSync('.').filter(f=>f.endsWith('.js')),...['anatomy','engine'].flatMap(dir=>readdirSync(dir).filter(f=>f.endsWith('.js')).map(f=>dir+'/'+f))];
const sourceHashes=Object.fromEntries(sourceFiles.map(f=>[f,createHash('sha256').update(readFileSync(f)).digest('hex')]));
writeFileSync(join(out,'source-manifest.json'),JSON.stringify({kind:'browser input files SHA256; source checkout must remain unchanged during audit',instrumentation:{exactSourceWeld:process.env.QA_SOURCE_MESH==='1'},files:sourceHashes},null,2)+'\n');
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
const server = process.env.QA_URL ? null : await staticServer();
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
if(shoulderPreview){
 await page.addInitScript(data=>{globalThis.__qaShoulderPreviewField=data},preview);
 await page.route("**/engine/shoulder-kinematics.js*",async route=>{
  const source=readFileSync("engine/shoulder-kinematics.js","utf8");
  const needle="if(!cached&&useField)a=root.ShoulderPoseField.sample(this.field,elev,plane);";
  assert(source.includes(needle));
  const patched=source.replace(needle,needle+`
  if(!cached&&useField&&elev>40&&elev<=180&&plane>=0&&plane<=90){
   const f=root.__qaShoulderPreviewField;
   if(!f||f.baseSignature!==this.fieldSignature)throw Error('Browser preview signature mismatch');
   const t=clamp((elev-40)/5,0,1),g=t*t*(3-2*t),d=root.ShoulderPoseField.sample(f,elev,plane);
   a=a.map((x,j)=>x+g*d[j]);
  }`).replace("backend:useField?'validated-field':'reference'","backend:useField?(elev>40&&elev<=180&&plane>=0&&plane<=90?'qa-additive-preview':'validated-field'):'reference'");
  await route.fulfill({body:patched,contentType:"text/javascript"});
 });
}
const instrument = String.raw`
let qaSkipDraw=false;const qaOriginalRender=render;render=()=>{if(!qaSkipDraw)qaOriginalRender();else{lastW=canvas.clientWidth;lastH=canvas.clientHeight;proj=ortho(state.zoom,lastW/lastH,.01,100);const eye=add(state.target,[18*Math.sin(state.yaw)*Math.cos(state.pitch),18*Math.sin(state.pitch),18*Math.cos(state.yaw)*Math.cos(state.pitch)]);view=lookAt(eye,state.target);pv=mm(proj,view)}};
const qaContactCurrent=v73Contact;let qaLegacyHead=false;
v73Contact=(spec,b)=>{const c=qaContactCurrent(spec,b);return qaLegacyHead&&c?.unilateral?{...c,unilateral:false}:c};
const qaExport=(f,id)=>{const B=f.binding,samples=[];for(let i=0;i<B.weights.length;i+=Math.max(1,Math.floor(B.weights.length/80))){const p=Array.from(B.rest.slice(i*3,i*3+3)),t=B.weights[i],l=B.lane[i],materialT=B.material[i],q=MuscleDeformation.point(f,p,t,l,materialT);samples.push({vertexId:i,p,t,l,materialT,q,originBlend:B.originBlend?.[i]||0})}return {id,name:B.definition.name,input:{side:v56Side,elev:settings.elev,plane:settings.plane,axial:settings.axial,elbow:settings.elbow,pronation:settings.pronation,sync:settings.sync},samples,state:{a:f.a,b:f.b,rest:Array.from(f.restFlat),now:Array.from(f.nowFlat),sheet:B.sheet,blendRange:B.blendRange,lanes:B.field.controls.length,tendonStart:B.definition.tissues.tendonStart,reference:B.reference,contact:f.contact,multiBone:!!f.originAt,originMix:f.originMix}}};
window.__qa={
 exportLevator(){
  const B=v58Meshes.get('v58Levator').common;
  return {rest:Array.from(B.rest),faces:Array.from(B.faces),weights:Array.from(B.weights),sourceOrigin:SOURCE_ORIGIN,scale:atlasScale,pins:B.definition.attachments.map(a=>Array.from(a.indices))};
 },
 levatorExperiment(){
  const spec=V60_ALL.find(s=>s.id==='v58Levator'),m=v58Meshes.get(spec.id),B=m.common;
  const attempt=()=>{try{return {candidate:window.__qaLevatorMasks?GeodesicAttachmentField.rebindSurface(B,window.__qaLevatorMasks):GeodesicAttachmentField.rebuild(B)}}catch(e){return {rejected:String(e)}}};
  if(!window.__qaLevatorCandidate)window.__qaLevatorCandidate=attempt();
  const {candidate,rejected}=window.__qaLevatorCandidate;
  function audit(side){
   const run=()=>{const f=v73MuscleState(spec,m,side),report=DeformationAudit.audit(f,MuscleDeformation.deform(f));return {report,trace:MuscleOutlierTrace.trace(f,report)}};
   const baseline=run();
   if(!candidate)return {side,baseline,rejected};
   try{m.common=candidate.binding;return {side,baseline,candidate:run(),bindingDiagnostics:candidate.diagnostics}}finally{m.common=B}
  }
  v56RefreshOther();
  return [audit(v56Side),v56WithOther(()=>audit(v56Inactive()))];
 },
 fullShoulderAudit(){
  v56RefreshOther();
  const side=()=>({side:v56Side,rig:window.__labShoulderV076(),muscles:V60_ALL.map(spec=>{const m=v58Meshes.get(spec.id),f=v73MuscleState(spec,m),B=f.binding;return {...DeformationAudit.audit(f,MuscleDeformation.deform(f)),vertices:B.weights.length,triangles:B.faces.length/3}})});
  // Active-side rig diagnostic flushes the pose; do not invoke it within
  // v56WithOther, whose snapshot is deliberately temporary.
  const active=side();
  const other=v56WithOther(()=>({side:v56Inactive(),muscles:V60_ALL.map(spec=>{const m=v58Meshes.get(spec.id),f=v73MuscleState(spec,m,v56Inactive()),B=f.binding;return {...DeformationAudit.audit(f,MuscleDeformation.deform(f)),vertices:B.weights.length,triangles:B.faces.length/3}})}));
  const neck=[...v73NeckMuscles.values()].map(m=>{const f=v73NeckMuscleState(m),B=f.binding;return {...DeformationAudit.audit(f,MuscleDeformation.deform(f)),vertices:B.weights.length,triangles:B.faces.length/3}});
  return {active,other,neck,glError:gl.getError(),arms:window.__labUpperBody().poses};
 },
 exportBindings(){return V60_ALL.map(spec=>{const m=v58Meshes.get(spec.id),B=m.common;return {spec,selection:m.selection,rest:Array.from(B.rest),faces:Array.from(B.faces),weights:Array.from(B.weights),field:{controls:B.field.controls},lane:Array.from(B.lane),material:Array.from(B.material),sheet:B.sheet,reference:B.reference,blendRange:B.blendRange,definition:{id:B.definition.id,name:B.definition.name,source:B.definition.source,tissues:B.definition.tissues,morphology:B.definition.morphology,attachmentFieldSource:B.definition.attachmentFieldSource}}})},
 rigFixture(){const model=v76ShoulderModel();return {landmarks:model.landmarks,contact:model.contact,scale:model.scale}},
 vertex(id,i){const spec=V60_ALL.find(s=>s.id===id),f=v73MuscleState(spec,v58Meshes.get(id)),B=f.binding;return {id,vertexId:i,point:MuscleDeformation.point(f,Array.from(B.rest.slice(i*3,i*3+3)),B.weights[i],B.lane[i],B.material[i]),a:f.a,b:f.b,weight:B.weights[i],lane:B.lane[i],materialT:B.material[i]}},
 ensureGPU(){if(!v73Renderer)v73Renderer=new MuscleRenderer(gl,shaderFS,(v,f)=>gpu(mesh(v,f)))},
 skipRender(value){qaSkipDraw=value;if(!value)request()},
 pose(p){for(const [k,v] of Object.entries(p)){$(k).step='any';$(k).value=v}update()},
 both(right,left){v56ChooseSide('right');this.pose(right);v56ChooseSide('left');this.pose(left);v56ChooseSide('right')},
 neck(n){v73SetNeck(n)},
 view(yaw){state.yaw=yaw;state.pitch=.1;fitPose();request()},
 legacyHead(value){qaLegacyHead=value;v56Other.key=null;v56RefreshOther();request()},
 camera(yaw,prior){if(prior)Object.assign(state,prior);else{state.yaw=yaw;state.pitch=.1;fitPose()}return {yaw:state.yaw,pitch:state.pitch,zoom:state.zoom,target:state.target.slice()}},
 headContactAudit(){return V60_ALL.filter(s=>MuscleDefinitions.configs[s.id]?.[1]==='head').map(s=>{const m=v58Meshes.get(s.id),f=v73MuscleState(s,m),B=f.binding,now=MuscleDeformation.deform(f),env=f.contact;let outward=0,inward=0,newPenetration=0,largestCorrection=0;for(let i=0;i<B.weights.length;i++){const p=Array.from(B.rest.slice(i*3,i*3+3)),q=Array.from(now.slice(i*3,i*3+3)),pre=MuscleDeformation.point({...f,contact:null},p,B.weights[i],B.lane[i],B.material[i]);const a=len(sub(pre,env.center)),b=len(sub(q,env.center)),rr=len(sub(p,env.restCenter)),c=(b-a)/atlasScale;if(c>1e-4)outward++;if(c< -1e-4)inward++;largestCorrection=Math.max(largestCorrection,Math.abs(c));newPenetration=Math.max(newPenetration,(Math.min(rr,env.radii[0])-b)/atlasScale)}return {id:s.id,outward,inward,largestCorrection,newPenetration,geometry:DeformationAudit.audit(f,now)}})},
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
  const originalHTML=readFileSync("index.html", "utf8");
  const sourceHTML=process.env.QA_SOURCE_MESH==="1"?originalHTML.replace("const key=p.map(x=>Math.round(x/grid)).join(',');","const key=p.join(',');"):originalHTML;
  const html = sourceHTML.replace(
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
  await page.goto(process.env.QA_URL || server.url);
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
  if(process.env.QA_EXPORT_LEVATOR==="1"){writeFileSync(join(out,"levator-binding.json"),JSON.stringify(await page.evaluate(()=>window.__qa.exportLevator())));await browser.close();if(server)await server.close();process.exit(0)}
  if(process.env.QA_EXPORT_BINDINGS==="1"){writeFileSync(join(out,"bindings.json"),JSON.stringify(await page.evaluate(()=>window.__qa.exportBindings()))); await browser.close(); if(server)await server.close(); process.exit(0)}
  await page.evaluate(() => window.__qa.ensureGPU());
  const initial = await page.evaluate(() => window.__qa.status());
  writeFileSync(join(out,"rig-calibration.json"),JSON.stringify(await page.evaluate(()=>window.__qa.rigFixture()),null,2));
  assert.equal(initial.backend, "gpu");
  assert.equal(initial.glError, 0);
  assert.deepEqual(
    (await page.evaluate(() => window.__qa.sample())).filter((x) => x.error),
    [],
    "missing common binding",
  );
  if(levatorExperiment){
    if(levatorMasks)await page.evaluate(m=>{window.__qaLevatorMasks=m},levatorMasks);
    await page.addScriptTag({content:readFileSync('scripts/experimental/geodesic-attachment-field.js','utf8')});
    await page.addScriptTag({content:readFileSync('scripts/experimental/trace-muscle-outliers.js','utf8')});
    await page.evaluate(()=>window.__qa.skipRender(true));
    const results=[];
    for(const [name,r,l,n] of [
      ['neutral',[0,0],[0,0],0],['low40',[40,40],[30,90],0],['recording70',[70,39],[0,0],0],
      ['side90',[90,0],[90,0],0],['high120',[120,40],[120,90],0],
      ['recording180',[180,44],[150,90],60],['asymmetric',[120,40,30,90,45],[75,90,-30,30,-45],-60],['return-neutral',[0,0],[0,0],0]
    ]){
      const pose=([elev,plane,axial=0,elbow=0,pronation=0])=>({elev,plane,axial,elbow,pronation});
      const sides=await page.evaluate(([r,l,n])=>{window.__qa.both(r,l);window.__qa.neck(n);return window.__qa.levatorExperiment()},[pose(r),pose(l),n]);
      results.push({name,right:pose(r),left:pose(l),neck:n,sides});
      writeFileSync(join(out,'levator-experiment.json'),JSON.stringify({results,preview:shoulderPreview,errors,productAccepted:false,anatomicalAcceptance:false},null,2));
      for(const s of sides){
       const compact=x=>x?{edge:x.report.edgeRatio,area:x.report.areaRatio,attachmentMm:x.report.attachmentMaxMm,degenerate:x.report.newDegenerateFaces,invalidFrames:x.report.invalidFrames,worstEdge:x.trace.worstEdge,worstArea:x.trace.worstArea,vertices:x.trace.vertices.map(v=>({i:v.i,weight:v.weight,lane:v.lane,materialT:v.materialT,restWidthMm:v.restWidthMm,restOffsetMm:v.restOffsetMm}))}:null;
       console.log('LEVATOR_TRACE '+JSON.stringify({name,side:s.side,baseline:compact(s.baseline),candidate:compact(s.candidate),binding:s.bindingDiagnostics,rejected:s.rejected,productAccepted:false}));
      }
    }
    assert.deepEqual(errors,[]);
  }else if(shoulderCompare){
    await page.evaluate(()=>window.__qa.skipRender(true));
    const cases=[
      {name:"neutral",right:{elev:0,plane:0},left:{elev:0,plane:0},neck:0},
      {name:"low40",right:{elev:40,plane:40},left:{elev:30,plane:90},neck:0},
      {name:"recording70",right:{elev:70,plane:39},left:{elev:0,plane:0},neck:0},
      {name:"side90",right:{elev:90,plane:0},left:{elev:90,plane:0},neck:0},
      {name:"high120",right:{elev:120,plane:40},left:{elev:120,plane:90},neck:0},
      {name:"recording180",right:{elev:180,plane:44},left:{elev:150,plane:90},neck:60},
      {name:"asymmetric",right:{elev:120,plane:40,axial:30,elbow:90,pronation:45},left:{elev:75,plane:90,axial:-30,elbow:30,pronation:-45},neck:-60},
      {name:"return-neutral",right:{elev:0,plane:0,axial:0,elbow:0,pronation:0},left:{elev:0,plane:0,axial:0,elbow:0,pronation:0},neck:0}
    ];
    const results=[];
    for(const input of cases){
      const audit=await page.evaluate(p=>{window.__qa.both(p.right,p.left);window.__qa.neck(p.neck);return window.__qa.fullShoulderAudit()},input);
      assert.equal(audit.active.muscles.length,17);assert.equal(audit.other.muscles.length,17);assert.equal(audit.neck.length,6);
      const expected=shoulderPreview&&input.right.elev>40?"qa-additive-preview":"validated-field";
      assert.equal(audit.active.rig.poseCanonicalRight.diagnostics.backend,expected);
      assert(audit.active.rig.connections.every(x=>x.sharedPivotErrorMm<1e-7));
      assert.equal(audit.glError,0);
      const all=[...audit.active.muscles,...audit.other.muscles,...audit.neck];
      assert(all.every(x=>x.nonfinite===0));
      results.push({input,...audit});
      writeFileSync(join(out,"shoulder-muscle-report.json"),JSON.stringify({preview:shoulderPreview,previewHash:shoulderPreview?createHash('sha256').update(readFileSync(previewFile)).digest('hex'):null,results,errors,anatomicalAcceptance:false,physicalDevice:false},null,2));
      console.log("SHOULDER_MUSCLE_CASE "+JSON.stringify({name:input.name,preview:shoulderPreview,backend:expected,maxEdge:Math.max(...all.map(x=>x.edgeRatio[1])),maxArea:Math.max(...all.map(x=>x.areaRatio[1])),maxAttachmentMm:Math.max(...all.map(x=>x.attachmentMaxMm??0)),degenerate:all.reduce((s,x)=>s+x.newDegenerateFaces,0),invalidFrames:all.reduce((s,x)=>s+x.invalidFrames,0)}));
    }
    assert.deepEqual(errors,[]);
    for(const [file,hash]of Object.entries(sourceHashes))assert.equal(createHash('sha256').update(readFileSync(file)).digest('hex'),hash,'Source changed during audit: '+file);
  }else{
  if(process.env.QA_HEAD_COMPARE==="1"){
    const result=[];
    await page.evaluate(()=>window.__qa.skipRender(false));
    for(const [elev,plane] of [[0,0],[90,0],[150,0],[180,0],[150,90],[180,90]]){
      await page.evaluate(([e,p])=>window.__qa.both({elev:e,plane:p},{elev:0,plane:0}),[elev,plane]);
      const cameras={};
      for(const legacy of [true,false]){
        await page.evaluate(x=>window.__qa.legacyHead(x),legacy);
        const metrics=await page.evaluate(()=>window.__qa.headContactAudit());
        result.push({elev,plane,legacy,metrics});
        for(const [name,yaw] of [["front",Math.PI],["side",Math.PI/2],["back",0]]){
          const cam=await page.evaluate(([y,c])=>window.__qa.camera(y,c),[yaw,cameras[name]||null]);
          if(cameras[name])assert.deepEqual(cam,cameras[name]);else cameras[name]=cam;
          await page.evaluate(()=>window.__qa.render());
          await page.screenshot({path:join(out,`${elev}-${plane}-${legacy?"before":"after"}-${name}.png`)});
        }
      }
    }
    writeFileSync(join(out,"head-contact-compare.json"),JSON.stringify(result,null,2));
    assert.equal(errors.length,0,errors.join("\n"));
    console.log("HEAD_CONTACT_COMPARE "+JSON.stringify(result.map(x=>({...x,metrics:x.metrics.map(({geometry,...m})=>m)}))));
  }
  await page.evaluate(()=>{window.__qa.legacyHead(false);window.__qa.skipRender(true)});
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
    maxReturnMm = 0,
    worstStep = null;
  for (const plane of [0, 40, 90]) {
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
            const p = prior[m].samples[i].q;
            const stepMm = Math.hypot(...q.map((x,k)=>x-p[k]))/.00555;
            if (stepMm > maxStepMm) {
              maxStepMm=stepMm;
              worstStep={stepMm,muscleId:samples[m].id,name:samples[m].name,vertexId:samples[m].samples[i].vertexId,plane,direction:"forward",from:{elev:elev-3,input:prior[m].input,sample:prior[m].samples[i],state:prior[m].state},to:{elev,input:samples[m].input,sample:samples[m].samples[i],state:samples[m].state}};
            }
          }
        }
      forward.push(samples);
      prior = samples;
    }
    for (let elev = 180; elev >= 0; elev -= 3) {
      await page.evaluate((p) => window.__qa.pose(p), { elev, plane });
      const samples = await page.evaluate(() => window.__qa.sample()),
        ref = forward[elev / 3];
      for (let m = 0; m < samples.length; m++)
        for (let i = 0; i < samples[m].samples.length; i++)
          maxReturnMm = Math.max(
            maxReturnMm,
            Math.hypot(
              ...samples[m].samples[i].q.map(
                (x, k) => x - ref[m].samples[i].q[k],
              ),
            ) / 0.00555,
          );
    }
  }
  writeFileSync(join(out,"worst-step.json"),JSON.stringify(worstStep,null,2));
  const refinement=[];
  for (const item of [worstStep,{muscleId:"v59TrapLower",vertexId:6696,plane:0,from:{elev:120},to:{elev:123}}]) {
    const rows=[];
    for(const step of [3,1.5,.75,.375]) {
      let last=null,maxMm=0;
      for(let elev=item.from.elev;elev<=item.to.elev+1e-8;elev+=step) {
        await page.evaluate(p=>window.__qa.pose(p),{elev,plane:item.plane});
        const v=await page.evaluate(([id,i])=>window.__qa.vertex(id,i),[item.muscleId,item.vertexId]);
        if(last)maxMm=Math.max(maxMm,Math.hypot(...v.point.map((x,k)=>x-last.point[k]))/.00555);
        last=v;
      }
      rows.push({stepDeg:step,maxMm});
    }
    const ratio=rows.at(-1).maxMm/Math.max(rows[0].maxMm,1e-9);
    refinement.push({id:item.muscleId,vertexId:item.vertexId,plane:item.plane,from:item.from.elev,to:item.to.elev,rows,ratio});
    assert(ratio<.3,"Step refinement found an unresolved discontinuity: "+JSON.stringify(refinement.at(-1)));
  }
  writeFileSync(join(out,"continuity-refinement.json"),JSON.stringify(refinement,null,2));
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
              contactUnilateral: c.unilateral ? 1 : 0,
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
  for(const [file,hash]of Object.entries(sourceHashes))assert.equal(createHash('sha256').update(readFileSync(file)).digest('hex'),hash,'Source changed during audit: '+file);
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
        worstStep: {id:worstStep.muscleId,name:worstStep.name,vertexId:worstStep.vertexId,plane:worstStep.plane,from:worstStep.from.elev,to:worstStep.to.elev},
        refinement,
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
  if (server) await server.close();
}
