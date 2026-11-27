// Reproduction of visible recording controls, not frame-by-frame tracking.
// Uses pinned full atlas in desktop Chromium emulating a narrow viewport.
// Not Safari or physical-device FPS certification.
import { chromium } from 'playwright';
import { readFileSync, writeFileSync, mkdirSync, existsSync } from 'node:fs';
import { join } from 'node:path';
import assert from 'node:assert/strict';
import { staticServer } from './helpers/static-server.mjs';
const out='test-results/v076-recording-replay',cache=process.env.ATLAS_CACHE||'/tmp/upper-limb-atlas';
mkdirSync(out,{recursive:true});mkdirSync(cache,{recursive:true});
const server=await staticServer();
const browser=await chromium.launch({headless:true,args:['--use-angle=swiftshader','--enable-unsafe-swiftshader','--ignore-certificate-errors']});
const page=await browser.newPage({viewport:{width:390,height:844},hasTouch:true,isMobile:true,deviceScaleFactor:1});
const errors=[],assetErrors=[],rows=[];
const pending=new Map();
page.on('pageerror',e=>errors.push(String(e)));
await page.route('**/*.stl',async route=>{
 const file=new URL(route.request().url()).pathname.split('/').at(-1),p=join(cache,file);
 try {
  if(!existsSync(p)){
   if(!pending.has(file))pending.set(file,(async()=>{
    const response=await fetch('https://raw.githubusercontent.com/Kevin-Mattheus-Moerman/BodyParts3D/2c83c698788cff57f6317bddcf712e7acfbea7d8/assets/BodyParts3D_data/stl/'+file);
    if(!response.ok)throw Error(response.status+' '+file);
    writeFileSync(p,Buffer.from(await response.arrayBuffer()));
   })());
   await pending.get(file);
  }
  await route.fulfill({body:readFileSync(p),contentType:'application/octet-stream',headers:{'access-control-allow-origin':'*'}});
 }catch(e){assetErrors.push(String(e));await route.abort();}
});
async function settle(){
 await page.waitForFunction(()=>!window.__labPerformance?.().pendingPose&&!window.__labPerformance?.().moving,null,{timeout:60000});
 await page.evaluate(()=>new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r))));
}
async function setPose(side,elev,plane){
 await page.evaluate(({side,elev,plane})=>{
  const select=document.getElementById('v56Side');select.value=side;select.dispatchEvent(new Event('change'));
  for(const [key,val] of Object.entries({elev,plane})){
   const el=document.getElementById(key);el.value=val;el.dispatchEvent(new Event('input',{bubbles:true}));
  }
  document.getElementById('elev').dispatchEvent(new Event('change',{bubbles:true}));
 },{side,elev,plane});
 await settle();
}
async function snapshot(name){
 const data=await page.evaluate(()=>{
  const camera=window.__labCamera(),bounds=window.__labBounds(),pose=window.__labUpperBody(),rig=window.__labShoulderV076();
  const [w,h]=camera.canvas,edge=document.body.classList.contains('edge-ui');
  const pads=edge?[Math.min(44,Math.max(20,w*.07)),Math.min(72,Math.max(52,w*.16)),Math.min(96,Math.max(64,h*.08)),Math.min(92,Math.max(56,h*.075))]:[Math.min(26,w*.065),Math.min(26,w*.065),Math.min(30,h*.065),Math.min(30,h*.065)];
  return {camera,bounds,pads,arms:pose.poses,rig:{backend:rig.poseCanonicalRight.diagnostics.backend,motion:rig.motion,sourceBoneSafeguard:rig.sourceBoneSafeguard,connections:rig.connections}};
 });
 const [L,R,T,B]=data.pads,s=data.bounds;
 assert(s.minX>=L-.2&&s.maxX<=s.width-R+.2&&s.minY>=T-.2&&s.maxY<=s.height-B+.2,JSON.stringify({name,...data}));
 rows.push({name,...data});
 await page.screenshot({path:join(out,name+'.png'),timeout:120000});
 console.log(JSON.stringify({name,bounds:s,pads:data.pads,zoom:data.camera.zoom,backend:data.rig.backend}));
}
try{
 await page.goto(server.url,{waitUntil:'domcontentloaded',timeout:120000});
 await page.waitForFunction(()=>window.__labDeepMuscles?.().loaded.length===17&&window.__labCervicalV073?.().neckMuscles===6,null,{timeout:240000});
 await page.evaluate(()=>document.getElementById('resetCam').click());await settle();
 assert((await page.evaluate(()=>window.__labCamera())).protectOnJointEnd);
 await snapshot('neutral');
 await setPose('right',70,39);await snapshot('right70-plane39');
 await setPose('right',180,44);await snapshot('right180-plane44');
 await page.evaluate(()=>{const el=document.getElementById('neckRotation');el.value=60;el.dispatchEvent(new Event('input',{bubbles:true}))});await settle();
 await snapshot('right180-neck60');
 await setPose('left',150,90);await snapshot('both-right180-left150');
 // Intentional magnification is not silently undone by inactivity.
 await page.evaluate(()=>document.getElementById('zoomIn').click());
 const manual=await page.evaluate(()=>window.__labCamera());assert(manual.manualFraming);
 await page.waitForTimeout(350);await settle();
 const after=await page.evaluate(()=>window.__labCamera());
 assert.equal(after.zoom,manual.zoom);assert.deepEqual(after.target,manual.target);
 await page.evaluate(()=>document.getElementById('resetCam').click());await settle();
 await snapshot('explicit-refit');
 assert.deepEqual(errors,[]);assert.deepEqual(assetErrors,[]);
}finally{
 writeFileSync(join(out,'report.json'),JSON.stringify({rows,errors,assetErrors,anatomicalAcceptance:false,physicalDevice:false,source:'visible recording joint values; neck60 is a test input, not recovered measurement'},null,2)+'\n');
 await browser.close();await server.close();
}
