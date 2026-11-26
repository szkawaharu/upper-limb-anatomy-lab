// All local browser modules are taken from the named git reference for the
// baseline. Each pair shares the exact camera and UI joint inputs. A changed
// calibration changes physical HT orientation; this is recorded, not hidden.
import { chromium } from 'playwright';
import { readFileSync,writeFileSync,mkdirSync,existsSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
import { join } from 'node:path';
import assert from 'node:assert/strict';
import { staticServer } from './helpers/static-server.mjs';
const reference='c81256f653c04c042ed917279b4fa3e4bde0f841',out=process.env.QA_OUTPUT||'test-results/comparison-v076',cache=process.env.ATLAS_CACHE||'/tmp/upper-limb-atlas';
mkdirSync(out,{recursive:true});
const server=await staticServer(),browser=await chromium.launch({executablePath:process.env.CHROME_PATH||'/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',headless:true,args:['--ignore-certificate-errors','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
const instrument=String.raw`
const comparisonRender=render;
render=()=>{lastW=canvas.clientWidth;lastH=canvas.clientHeight;proj=ortho(state.zoom,lastW/lastH,.01,100);const eye=add(state.target,[18*Math.sin(state.yaw)*Math.cos(state.pitch),18*Math.sin(state.pitch),18*Math.cos(state.yaw)*Math.cos(state.pitch)]);view=lookAt(eye,state.target);pv=mm(proj,view)};
const comparisonVisibility=new Map([...new Set([...muscleToggleIds,'v73NeckMusclesVisible','v56Back'])].map(id=>[id,$(id)?.checked]));
window.__compare={
 set(p,mode,cam,yaw){
  for(const [id]of comparisonVisibility)if($(id))$(id).checked=mode!=='bones';
  for(const side of ['right','left']){v56ChooseSide(side);for(const k of ['elev','plane','axial','elbow','pronation','wrist','grip'])if($(k)){$(k).value=p[side]?.[k]||0}update()}
  v56ChooseSide('right');v73SetNeck(p.neck||0);
  if(cam)Object.assign(state,cam);else{state.yaw=yaw;state.pitch=.1;fitPose()}
  request();return {yaw:state.yaw,pitch:state.pitch,zoom:state.zoom,target:state.target.slice()};
 },
 draw(){comparisonRender()},
 audit(){return {rig:window.__labShoulderV076?.(),arms:window.__labUpperBody().poses,neck:window.__labCervicalV073(),backend:v73Renderer?.backend,glError:gl.getError()}},
 diagnostics(){for(const id of ['rigAxes','rigLandmarks','rigSurface','rigAttachments'])if($(id)){$(id).checked=true;$(id).dispatchEvent(new Event('change'))}request();comparisonRender()}
};`;
const modes=(process.env.QA_COMPARE_MODES||'bones,muscles').split(',');
const previous=process.env.QA_COMPARE_RESUME==='1'?JSON.parse(readFileSync(join(out,'comparison.json'))):[];
const results=previous.filter(r=>!modes.includes(r.mode)),cameras=new Map(results.map(r=>[[r.name,r.mode,r.view].join('-'),r.camera]));
const poses=[['neutral',{}],['side90',{right:{elev:90}}],['front150',{right:{elev:150,plane:90}}],['asymmetric',{right:{elev:120,plane:40,axial:30,elbow:90,pronation:45},left:{elev:75,plane:90,axial:-30,elbow:30,pronation:-45},neck:60}]];
try{
 for(const version of ['v075','v076']){
  const page=await browser.newPage({viewport:{width:1280,height:900}}),errors=[];
  page.on('pageerror',e=>{errors.push(String(e));console.log('PAGE ERROR',version,String(e))});
  await page.route('**/*.stl',async route=>{
   const file=new URL(route.request().url()).pathname.split('/').at(-1),p=join(cache,file);
   assert(existsSync(p),'Missing pinned atlas cache '+file);
   await route.fulfill({body:readFileSync(p),contentType:'application/octet-stream',headers:{'access-control-allow-origin':'*'}});
  });
  if(version==='v075')await page.route(/\/[^?]+\.js(?:\?.*)?$/,async route=>{
   const p=decodeURIComponent(new URL(route.request().url()).pathname).slice(1);
   await route.fulfill({body:execFileSync('git',['show',reference+':'+p]),contentType:'text/javascript'});
  });
  const source=version==='v075'?execFileSync('git',['show',reference+':index.html'],{encoding:'utf8'}):readFileSync('index.html','utf8');
  await page.route('**/index.html',route=>route.fulfill({body:source.replace(/\}\)\(\);\s*<\/script>/,()=>instrument+'})();</script>'),contentType:'text/html'}));
  await page.goto(server.url,{waitUntil:'domcontentloaded',timeout:120000});
  await page.waitForFunction(()=>window.__labDeepMuscles?.().loaded.length===17&&window.__labCervicalV073?.().neckMuscles===6,null,{timeout:240000});
  console.log('Atlas loaded',version);
  for(const [name,pose]of poses)for(const mode of modes)for(const [view,yaw]of [['front',Math.PI],['side',Math.PI/2],['back',0],['oblique',Math.PI*.75]]){
   const key=[name,mode,view].join('-'),camera=await page.evaluate(([p,m,c,y])=>window.__compare.set(p,m,c,y),[pose,mode,cameras.get(key),yaw]);
   if(!cameras.has(key))cameras.set(key,camera);else assert.deepEqual(camera,cameras.get(key));
   await page.evaluate(()=>window.__compare.draw());
   const audit=await page.evaluate(()=>window.__compare.audit());
   assert.equal(audit.glError,0);assert.deepEqual(errors,[]);
   if(version==='v076')assert.equal(audit.rig.poseCanonicalRight.diagnostics.backend,'validated-field','Unexpected CPU rig fallback: '+JSON.stringify(audit.rig));
   await page.screenshot({clip:await page.locator('#viewer').boundingBox(),path:join(out,version+'-'+key+'.png'),timeout:120000});
   results.push({version,reference:version==='v075'?reference:'working-tree',name,mode,visibility:mode==='bones'?'all muscle toggles off':'all muscle toggles including v56Back and neck on',view,pose,camera,audit,errors});
   writeFileSync(join(out,'comparison.json'),JSON.stringify(results,null,2));console.log(version,key);
  }
  if(version==='v076'){
   await page.evaluate(()=>window.__compare.diagnostics());
   await page.screenshot({clip:await page.locator('#viewer').boundingBox(),path:join(out,'v076-diagnostics.png'),timeout:120000});
  }
  await page.close();
 }
}finally{await browser.close();await server.close()}
