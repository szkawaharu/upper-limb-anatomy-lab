// Test-only static continuum snapshot vs identical source input to the curve
// engine. Explicit target pose; neither ROM nor production GPU is certified.
import{chromium}from'playwright';
import{readFileSync,writeFileSync,mkdirSync,existsSync}from'node:fs';
import{gunzipSync}from'node:zlib';
import{join}from'node:path';
import{createHash}from'node:crypto';
import assert from'node:assert/strict';
import{staticServer}from'./helpers/static-server.mjs';
const packetFile=process.env.QA_SOLID_INPUT||'test-results/v076-volumetric/neutral-preview.json.gz',packet=JSON.parse(gunzipSync(readFileSync(packetFile))),out=process.env.QA_OUTPUT||'test-results/v076-solid-neutral',cache=process.env.ATLAS_CACHE||'/tmp/upper-limb-atlas';
const cameraReference=process.env.QA_CAMERA_REFERENCE?JSON.parse(readFileSync(process.env.QA_CAMERA_REFERENCE,'utf8')):null;
assert.equal(packet.accepted,false);assert.equal(packet.report.converged,true);assert(Number.isFinite(packet.pose.elev)&&Number.isFinite(packet.pose.plane));mkdirSync(out,{recursive:true});
const instrument=String.raw`
const solidDraw=v73DrawMuscle,solidVisible=v58Visible,solidRender=render;let solidMode=null,solidIsolation=false,solidGeometries=null,solidAudit=null,solidId=null;
v58Visible=(spec)=>solidMode&&solidIsolation?!v56DrawMirror&&spec.id===solidId:solidVisible(spec);
// The offline snapshot is the RIGHT side only. Never reflect that snapshot
// onto the independently posed left arm (neutral happened to conceal this).
v73DrawMuscle=(spec,m)=>{if(spec.id!==solidId||!solidMode||v56DrawMirror)return solidDraw(spec,m);gl.useProgram(program);drawPart(solidGeometries[solidMode],ID,OR,[.67,.26,.30],1,1);return true;};
// Freeze scheduling while retaining the actual renderer for each capture.
render=()=>{lastW=canvas.clientWidth;lastH=canvas.clientHeight;proj=ortho(state.zoom,lastW/lastH,.01,100);const eye=add(state.target,[18*Math.sin(state.yaw)*Math.cos(state.pitch),18*Math.sin(state.pitch),18*Math.cos(state.yaw)*Math.cos(state.pitch)]);view=lookAt(eye,state.target);pv=mm(proj,view)};
window.__solid={
 install(packet){
  solidId=packet.id;
  for(const side of ['right','left']){v56ChooseSide(side);for(const key of ['elev','plane','axial','elbow','pronation'])$(key).value=0;update()};v56ChooseSide('right');v73SetNeck(0);for(const key of ['elev','plane'])$(key).value=packet.pose[key];update();
  const raw=packet.binding,spec=V60_ALL.find(x=>x.id===packet.id),p=Array.from({length:raw.rest.length/3},(_,i)=>raw.rest.slice(i*3,i*3+3)),a=v58Anchor(spec.a),b=v58Anchor(spec.b),definition=MuscleDefinitions.create(spec,p,raw.faces,raw.weights,[a,b]),B={...raw,definition,rest:Float32Array.from(raw.rest),faces:Uint32Array.from(raw.faces),weights:Float32Array.from(raw.weights),lane:Float32Array.from(raw.lane),material:Float32Array.from(raw.material),cache:new Map()},f=MuscleDeformation.pose(B,a,b,v73Contact(spec,b)),now=MuscleDeformation.deform(f),geometry=positions=>gpu(mesh(Array.from({length:positions.length/3},(_,i)=>Array.from(positions.slice(i*3,i*3+3))),Array.from(B.faces)));
  let traceErrorMm=null;if(packet.traceExpected){const t=packet.traceExpected;traceErrorMm=Math.hypot(...t.point.map((v,k)=>v-now[3*t.vertex+k]))/atlasScale;}
  const S=packet.surface||raw,C={rest:Float64Array.from(S.rest),faces:Uint32Array.from(S.faces),weights:Float64Array.from(S.weights),definition:B.definition};
  solidGeometries={curve:geometry(now),continuum:gpu(mesh(Array.from({length:packet.positions.length/3},(_,i)=>packet.positions.slice(i*3,i*3+3)),Array.from(C.faces)))};
  let error=0;for(let i=0;i<C.weights.length;i++)if(C.weights[i]===0||C.weights[i]===1){const target=AnatomyMath.apply(C.weights[i]===0?a:b,C.rest.slice(3*i,3*i+3));error=Math.max(error,Math.hypot(...target.map((v,k)=>v-packet.positions[3*i+k]))/atlasScale)};
  const volume=(positions,faces)=>DeformationAudit.signedVolume(positions,faces);
  solidAudit={a,b,traceErrorMm,attachmentMaxMm:error,sourceVertices:B.weights.length,sourceTriangles:B.faces.length/3,displayVertices:C.weights.length,displayTriangles:C.faces.length/3,curve:DeformationAudit.audit(f,now),continuum:DeformationAudit.audit({...f,binding:C,vertexDiagnostics:()=>({})},Float64Array.from(packet.positions)),float32:DeformationAudit.audit({...f,binding:{...C,rest:Float32Array.from(C.rest),topology:null},vertexDiagnostics:()=>({})},Float32Array.from(packet.positions)),curveSignedVolumeRatio:volume(now,B.faces)/volume(B.rest,B.faces),continuumSignedVolumeRatio:volume(packet.positions,C.faces)/volume(C.rest,C.faces)};return solidAudit;
 },
 set(mode,isolation,yaw,camera){solidMode=mode;solidIsolation=isolation;for(const id of muscleToggleIds)if($(id))$(id).checked=!isolation;if($('v56Back'))$('v56Back').checked=!isolation;update();if(camera)Object.assign(state,camera);else{state.yaw=yaw;state.pitch=.1;fitPose();state.zoom*=1.05}solidRender();return{yaw:state.yaw,pitch:state.pitch,zoom:state.zoom,target:state.target.slice()};},
 audit(){return{...solidAudit,rig:window.__labShoulderV076(),glError:gl.getError()}}
};`;
const server=await staticServer(),browser=await chromium.launch({executablePath:process.env.CHROME_PATH||'/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',headless:true,args:['--ignore-certificate-errors','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
try{
 const page=await browser.newPage({viewport:{width:1280,height:900}}),errors=[];page.on('pageerror',e=>errors.push(String(e)));
 await page.route('**/*.stl',route=>{const file=join(cache,new URL(route.request().url()).pathname.split('/').at(-1));assert(existsSync(file));return route.fulfill({body:readFileSync(file),contentType:'application/octet-stream',headers:{'access-control-allow-origin':'*'}})});
 const html=readFileSync('index.html','utf8');await page.route('**/index.html',route=>route.fulfill({body:html.replace(/\}\)\(\);\s*<\/script>/,()=>instrument+'})();</script>'),contentType:'text/html'}));
 await page.goto(server.url,{waitUntil:'domcontentloaded',timeout:120000});await page.waitForFunction(()=>window.__labDeepMuscles?.().loaded.length===17&&window.__labCervicalV073?.().neckMuscles===6,null,{timeout:240000});
 console.log('Atlas and cervical meshes loaded');const initial=await page.evaluate(p=>window.__solid.install(p),packet);assert(initial.attachmentMaxMm<1e-5,JSON.stringify(initial));assert(Number.isFinite(initial.curveSignedVolumeRatio)&&Number.isFinite(initial.continuumSignedVolumeRatio),'Invalid volume measurement');assert(Math.abs(initial.continuumSignedVolumeRatio-packet.report.volumeRatio)<1e-5,'Independent solid and browser surface volumes differ');
 for(const bone of ['a','b'])for(const key of ['q','p','rest'])for(let k=0;k<packet.target[bone][key].length;k++)assert(Math.abs(initial[bone][key][k]-packet.target[bone][key][k])<1e-9,'Final bone pose differs '+bone+key+k);
 if(packet.traceExpected)assert(initial.traceErrorMm<.001,'Captured production trace was not reproduced: '+initial.traceErrorMm+'mm');
 const rows=[];const initialAudit=await page.evaluate(()=>window.__solid.audit());assert.equal(initialAudit.glError,0);assert.deepEqual(errors,[]);
 if(!process.env.QA_NUMERIC_ONLY)for(const isolation of process.env.QA_ISOLATED_ONLY?[true]:[false,true])for(const [view,yaw]of [['front',Math.PI],['side',Math.PI/2],['back',0]]){
  let camera=cameraReference?.rows.find(row=>row.isolation===isolation&&row.view===view)?.camera||null;
  if(cameraReference)assert(camera,'Missing reference camera '+view);
  for(const mode of ['curve','continuum']){
   const actual=await page.evaluate(args=>window.__solid.set(...args),[mode,isolation,yaw,camera]);if(camera)assert.deepEqual(actual,camera);else camera=actual;
   const audit=await page.evaluate(()=>window.__solid.audit());assert.equal(audit.glError,0);assert.deepEqual(errors,[]);
   const file=mode+'-'+(isolation?'isolated':'all')+'-'+view+'.png';await page.screenshot({clip:await page.locator('#viewer').boundingBox(),path:join(out,file),timeout:120000});rows.push({mode,isolation,view,camera,file,audit});console.log(file);
  }
 }
 writeFileSync(join(out,process.env.QA_NUMERIC_ONLY?'audit-only.json':'comparison.json'),JSON.stringify({experimental:true,accepted:false,muscleId:packet.id,restSurfaceIdentical:packet.restSurfaceIdentical??true,pose:packet.pose,kind:packet.kind||'full-fem',boundaryComparison:packet.boundaryComparison||{identical:true},packetSha256:createHash('sha256').update(readFileSync(packetFile)).digest('hex'),inputs:packet.inputs,limitations:['static target only',packet.restSurfaceIdentical===false?'same source atlas with different mesh preparations; repaired source geometry is separately audited':'same repaired rest surface; continuum may use an exact surface subdivision','footprint comparison declared separately; all are estimates','no collision or smartphone performance certification','other muscles still production curve engine'],initialAudit,rows,errors},null,2)+'\n');
}finally{await browser.close();await server.close()}
