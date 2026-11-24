// v0.56 integration regression; offline synthetic/Canvas2D harness, not clinical or browser validation.
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
const html=readFileSync('index.html','utf8');
let js=html.match(/<script>([\s\S]*?)<\/script>/)[1];
new Function(js);
assert.ok(/上肢3D解剖ラボ · v\d+\.\d+/.test(html));
const elements=new Map(),errors=[],callbacks=[],generic=new Map();
function el(id=""){
 const attrs={},classes=new Set(),events={};
 return {id,value:"",checked:false,min:"0",max:"180",clientWidth:390,clientHeight:740,width:390,height:740,offsetWidth:56,style:{setProperty(){}},dataset:{},children:[],
 classList:{add(...a){a.forEach(x=>classes.add(x))},remove(...a){a.forEach(x=>classes.delete(x))},contains(x){return classes.has(x)},toggle(x,b){if(b===undefined)b=!classes.has(x);b?classes.add(x):classes.delete(x);return b}},
 addEventListener(n,f){(events[n]??=[]).push(f)},dispatch(n){for(const f of events[n]||[])f({target:this,preventDefault(){}})},click(){this.dispatch("click")},
 getAttribute(n){return attrs[n]??null},setAttribute(n,v){attrs[n]=String(v)},removeAttribute(n){delete attrs[n]},appendChild(x){this.children.push(x);return x},insertBefore(x){this.children.unshift(x);return x},
 append(...x){this.children.push(...x)},replaceChildren(...x){this.children=x},querySelector(s){return s.includes("viewer-toolbar")?get(s):null},querySelectorAll(){return[]},focus(){},scrollIntoView(){},
 getBoundingClientRect(){return{width:390,height:740,left:0,top:0}},getContext(t){if(t!=="2d")return null;return new Proxy({createLinearGradient(){return{addColorStop(){}}}},{get(o,k){return k in o?o[k]:()=>{}},set(o,k,v){o[k]=v;return true}})}}
}
const get=s=>{if(!generic.has(s))generic.set(s,el(s));return generic.get(s)};
for(const m of html.matchAll(/<([a-z][\w-]*)\b([^>]*\bid="([^"]+)"[^>]*)>/gi)){
 const [_,tag,attr,id]=m,e=el(id);
 for(const a of attr.matchAll(/([\w-]+)="([^"]*)"/g)){e.setAttribute(a[1],a[2]);if(["value","min","max"].includes(a[1]))e[a[1]]=a[2]}
 e.checked=/\bchecked(?:\s|$)/.test(attr);
 if(tag==="select"){const rest=html.slice(m.index+m[0].length),sel=rest.slice(0,rest.indexOf("</select>"));e.value=sel.match(/<option\b[^>]*value="([^"]+)"/)?.[1]||"";const selected=sel.match(/<option\b[^>]*value="([^"]+)"[^>]*selected/);if(selected)e.value=selected[1]}
 elements.set(id,e);
}
const doc={getElementById(id){return elements.get(id)||null},body:el("body"),querySelector:get,querySelectorAll(){return[]},createElement:el,createElementNS:(_,t)=>el(t),addEventListener(){},activeElement:el()};
const win={addEventListener(){},devicePixelRatio:1,scrollTo(){}};
const instrument=String.raw`
window.__runV56Tests=async()=>{
 const passed=[],ok=(b,m)=>{if(!b)throw Error(m);passed.push(m)};
 ok(Boolean(state.scene),'startup fallback scene');
 const data=window.__labUpperBody();
 ok(data.ribPairs===12&&v56Torso().parts.filter(p=>p.id.startsWith('rib-')).length===24,'24 rib meshes');
 ok(v56Torso().parts.filter(p=>p.id.startsWith('vertebra-')).length===19,'7 cervical + 12 thoracic vertebrae');
 ok(geomDynamic!==v56Other.geom,'separate dynamic buffer dictionaries');
 ok(len(sub(v56Reflect(v56Reflect([1,2,3])),[1,2,3]))<1e-9,'reflection involution');
 const a=transform(qa([0,1,0],.6),[.3,.7,1],[.2,-.4,.1]);
 const mat=mm(v56ReflectionMatrix(),mat4(qa([0,1,0],.6),[.3,.7,1]));
 const p=[.2,-.4,.1,1],m=[0,1,2].map(r=>p.reduce((sum,v,k)=>sum+v*mat[k*4+r],0));
 ok(len(sub(m,v56Reflect(a)))<1e-9,'WebGL mirrored model matrix');
 $('elbow').value='75';$('elev').value='35';update();v56ChooseSide('left');
 ok(settings.elbow===0&&settings.elev===0,'left remains neutral after right edit');
 $('elbow').value='110';$('pronation').value='35';$('wristFlex').value='20';update();v56ChooseSide('right');
 ok(settings.elbow===75&&settings.elev===35&&settings.pronation===0,'right pose restored');
 const s=window.__labUpperBody();
 ok(s.poses.left.elbow===110&&s.poses.left.pronation===35&&s.poses.left.wristFlex===20,'left independent elbow, rotation and wrist');
 const old=v56Other.key;$('elbow').value='90';update();v56RefreshOther();
 ok(v56Other.key===old,'inactive arm not rebuilt on active angle change');
 const pts=poseBoundingPoints();ok(pts.length>0&&pts.every(p=>p.every(Number.isFinite)),'finite bilateral camera bounds');
 fitPose();ok(state.zoom>0&&Number.isFinite(state.zoom),'camera fits both arms');
 render();ok(drawCalls.some(c=>c.mirror)&&drawCalls.some(c=>!c.mirror),'software rendering submits mirrored and original geometry');
 $('v56Both').checked=false;render();ok(!drawCalls.some(c=>c.mirror),'single arm hides mirror');
 v56ChooseSide('left');render();ok(drawCalls.some(c=>c.mirror),'single left arm remains mirrored');
 $('v56Both').checked=true;$('v56Copy').click();
 ok(window.__labUpperBody().poses.right.elbow===110,'copy pose to other arm');
 $('v56Back').checked=false;render();ok(!v56BackVisible(),'back toggle off');
 $('rib').checked=false;update();ok(!window.__labUpperBody().thorax,'thorax toggle off');
 setDisplayPreset('skeleton');ok(!$('v56Back').checked&&muscleToggleIds.every(id=>!$(id).checked),'skeleton preset hides back and arm muscles');
 $('restoreDisplayState').click();ok(!$('v56Back').checked,'display restore preserves previous back state');
 // Synthetic atlas-local meshes exercise the real-ready pose path without network.
 const synthetic={scapula:scapula(),clavicle:tube([[0,0,0,.06],[.65,.04,-.02,.06]]),
 humerus:tube([[0,0,0,.16],[.12,-.40,.02,.12],[.14,-1.45,0,.13]],{sides:16,steps:12}),
 radius:tube([[.13,0,0,.08],[.11,-.6,.03,.055],[.10,-1.3,.04,.10]]),
 ulna:tube([[-.13,0,0,.08],[-.13,-.6,0,.06],[-.12,-1.3,0,.065]])};
 for(const [id,g] of Object.entries(synthetic)){
  realSurfaceVertices[id]=g.vertices;realSurfaceFaces[id]=g.faces;realMesh[id]=gpu(g);
  realMesh[id].sourceVertices=g.vertices;realMesh[id].sourceFaces=g.faces;
 }
 resetScapularSurfaceProjection();realReady=true;setDisplayPreset('on');
 $('rib').checked=true;update();v56RefreshOther();v56BuildBack();
 ok(realPose.humerus&&v56Other.real.humerus,'synthetic atlas-ready pose for both arms');
 ok(v56BackIds.every(id=>geomDynamic[id]?.cpu?.vertices.every(p=>p.every(Number.isFinite))),'finite real-ready back meshes');
 ok(poseBoundingPoints().every(p=>p.every(Number.isFinite)),'finite real-ready bilateral bounds');
 ok(geomDynamic!==v56Other.geom&&geomDynamic.bicepsLong!==v56Other.geom.bicepsLong,'independent real-ready muscle buffers');
 const arm=realPose.humerus[1].slice();v56RefreshOther();
 ok(len(sub(realPose.humerus[1],arm))<1e-9,'inactive refresh preserves active joint center');

 ok(V57_MANIFEST.length===79&&new Set(V57_MANIFEST.map(p=>p.file)).size===79,'79 unique atlas files');
 ok(V57_MANIFEST.filter(p=>p.group==='ribs'&&p.side==='right').length===12&&V57_MANIFEST.filter(p=>p.group==='ribs'&&p.side==='left').length===12,'24 independently sourced ribs');
 ok(V57_MANIFEST.filter(p=>p.group==='spine').length===19,'19 individually sourced vertebrae');
 ok(V57_MANIFEST.filter(p=>p.group==='sternum').length===3&&V57_MANIFEST.filter(p=>p.group==='cartilage').length===16,'sternum and cartilage manifest');
 const fixture=new ArrayBuffer(134),dv=new DataView(fixture);dv.setUint32(80,1,true);
 const triangle=[[-14,-72,1341],[-24,-72,1341],[-14,-82,1341]];
 triangle.forEach((v,i)=>v.forEach((n,j)=>dv.setFloat32(96+i*12+j*4,n,true)));
 const parsed=v57ParseSTL(fixture,{id:'synthetic',bytes:134});
 ok(parsed.faces.length===3&&len(sub(parsed.vertices[0],atlasToWorld(triangle[0])))<1e-10,'binary STL common atlas coordinates');
 ok(len(sub(parsed.vertices[1],atlasToWorld(triangle[1])))<1e-10,'STL preserves absolute position and scale');
 let rejected=false;try{v57ParseSTL(fixture,{id:'bad',bytes:135})}catch(e){rejected=true}
 ok(rejected,'mismatched STL size rejected');
 const stage=V57_MANIFEST.map(e=>({...e,g:gpu(parsed),color:cBone}));
 ok(v57ValidateSet(stage)&&!v57ValidateSet(stage.slice(1)),'atomic complete set required');
 v57State.parts=stage;v57State.mode='atlas';v57State.loaded=79;
 $('v57TorsoLayer').value='ribs';ok(v57VisibleTorso().length===24,'atlas ribs filter');
 $('v57TorsoLayer').value='spine';ok(v57VisibleTorso().length===19,'atlas spine filter');
 $('v57TorsoLayer').value='all';ok(v57VisibleTorso().length===79,'complete atlas rendering set');
 ok(poseBoundingPoints().every(p=>p.every(Number.isFinite)),'atlas camera bounds finite');
 v57State.mode='fallback';v57State.parts=[];$('v57TorsoLayer').value='ribs';
 ok(v57VisibleTorso().length===24&&v57VisibleTorso().every(p=>p.id.startsWith('rib-')),'fallback group selection');
 $('v57TorsoLayer').value='all';


 ok(V571_CORE.length===48&&V571_COSTAL.length===14&&V571_DISCS.length===17,'core and augmentation groups');
 ok(V571_COSTAL.every(p=>p.level>=1&&p.level<=7),'costal levels 1–7 only');
 ok(!V571_DISCS.some(p=>p.id==='discC1'||p.id==='discT12'),'no C1–C2 or T12–L1 disc');
 ok(V571_COSTAL_COLOR!==V571_DISC_COLOR,'cartilage and disc color differentiation');
 const make=e=>({...e,g:gpu(parsed),color:v571Color(e)});
 // Failure in new cartilage must not remove the precise base or discs.
 v57State.parts=[];v57State.mode='fallback';
 await v57LoadTorso(async e=>{if(e.id==='cartR3')throw Error('simulated missing cartilage');return make(e)});
 ok(v57State.mode==='partial'&&v571HasCore()&&v57State.parts.length===65,'optional cartilage failure retains core and discs');
 let retried=0;
 await v57LoadTorso(async e=>{retried++;return make(e)});
 ok(retried===14&&v57State.mode==='atlas'&&v57State.parts.length===79,'retry only incomplete group');
 $('v57TorsoLayer').value='discs';ok(v57VisibleTorso().length===17,'disc filter');
 $('v57TorsoLayer').value='cartilage';ok(v57VisibleTorso().length===16,'cartilage filter');
 $('v571LegacyCart').checked=false;ok(v57VisibleTorso().length===14,'old cartilage independent visibility');
 $('v571LegacyCart').checked=true;$('v57TorsoLayer').value='all';
 ok(v571BoxHit([[0,0,0],[1,1,1]],[[.5,.5,.5],[2,2,2]])&&!v571BoxHit([[0,0,0],[1,1,1]],[[2,2,2],[3,3,3]]),'cartilage overlap broad-phase');
 const j=v571JointAudit();ok(j.length===2&&j.every(x=>x.scToClavicleMm!==null&&x.acToScapulaMm!==null),'bilateral inferred shoulder-center audit');
 ok(v571MinDistance([0,0,0],[])===null,'missing geometry not reported as zero clearance');
 v57State.parts=[];v57State.mode='fallback';
 await v57LoadTorso(async e=>{if(e.id==='R1')throw Error('simulated missing base');return make(e)});
 ok(v57State.mode==='fallback'&&v57State.parts.length===0,'base failure keeps old fallback');

 {
 ok(v572Scale('auto',true,2)===.85&&v572Scale('auto',false,2)===2,'adaptive resolution');
 ok(v572Scale('high',true,2)===2&&v572Scale('light',false,2)===1,'manual quality overrides');
 ok(v572Scale('auto',true,.7)===.7,'adaptive resolution never exceeds device');
 v572Preset('forward');
 let v572p=window.__labUpperBody().poses;
 ok(v572p.right.elev===90&&v572p.left.elev===90&&v572p.right.plane===90&&v572p.left.plane===90,'bilateral forward preset');
 ok(v572p.right.elbow===0&&v572p.left.elbow===0&&v572p.right.wristFlex===0&&v572p.left.fingerCurl===0,'forward preset straight elbows and hands');
 ok(qr(shoulderPoseFor(90,90,0,true).qH,[0,-1,0])[2]<-.999999,'real anterior direction');
 const palm=v572ChoosePalm();ok(Number.isFinite(palm.alignment)&&palm.pronation>=-90&&palm.pronation<=90,'bounded palm optimization');
 v572Preset('t');v572p=window.__labUpperBody().poses;ok(v572p.right.plane===0&&v572p.left.elev===90,'bilateral T preset');
 v572Preset('up');v572p=window.__labUpperBody().poses;ok(v572p.right.elev===180&&v572p.left.elev===180,'bilateral overhead preset');
 v572Preset('neutral');v572p=window.__labUpperBody().poses;ok(v572p.right.elev===0&&v572p.left.elev===0,'bilateral neutral preset');
 $('elbow').value='42';v572QueuePose();v572QueuePose();
 ok(v572Pending&&settings.elbow===0,'pose changes coalesced before frame');
 v572FinishInteraction();
 ok(!v572Pending&&!v572Moving&&settings.elbow===42,'final pose flushed and quality restored');
 state.zoom=.75;$('v572Fit').checked=true;
 ok(v572FitIfNeeded()&&state.zoom>.75,'clipped model fit');
 $('v572Fit').checked=false;state.zoom=.75;
 ok(!v572FitIfNeeded()&&state.zoom===.75,'opt-out preserves manual zoom');
 v572Preset('neutral');
 }
 return passed;
};`;
js=js.replace(/\}\)\(\);\s*$/,instrument+"\n})();");
new Function("document","window","fetch","requestAnimationFrame","matchMedia","innerWidth","innerHeight","setTimeout","clearTimeout","console",js)(doc,win,async()=>{throw Error("offline test")},f=>callbacks.push(f),()=>({matches:true}),390,740,()=>0,()=>{}, {log(){},warn(){},error(...a){errors.push(a.map(String).join(" "))}});
assert.equal(errors.length,0,errors.join("\n"));
const results=await win.__runV56Tests();
console.log('PASS v0.57.2: '+results.length+' integration checks — '+results.join('; '));
