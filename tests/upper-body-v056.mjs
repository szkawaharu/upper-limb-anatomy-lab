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

 {
 ok(V58_DEEP.length===5&&new Set(V58_DEEP.map(x=>x.file)).size===5,'five distinct source muscles');
 ok(v56BackIds.length===4&&!v56BackIds.includes('v56Serratus'),'old deep plates removed');
 const box=[[0,0,0],[1,1,1]],spec=V58_DEEP[0];
 ok(v58Weight(spec,[0,.5,.5],box)===0&&v58Weight(spec,[1,.5,.5],box)===1,'rhomboid endpoint weights');
 ok(v58Weight(V58_DEEP[2],[.5,1,.5],box)===0&&v58Weight(V58_DEEP[2],[.5,0,.5],box)===1,'levator downward attachment weights');
 ok(v58Weight(V58_DEEP[3],[.5,.5,0],box)===0&&v58Weight(V58_DEEP[3],[.5,.5,1],box)===1,'serratus posterior attachment proxy');
 const a={q:ID,p:OR,rest:OR},b={q:qa([0,1,0],.4),p:[1,2,3],rest:[.1,.2,.3]},point=[.4,.5,.6];
 ok(len(sub(v58Blend(a,b,point,0),point))<1e-10&&len(sub(v58Blend(a,b,point,1),v58Transform(b,point)))<1e-10,'skin endpoints follow respective anchor');
 const matrix=v58AnchorMatrix(b),mp=[...point,1],mv=[0,1,2].map(r=>mp.reduce((s,v,k)=>s+v*matrix[k*4+r],0));
 ok(len(sub(mv,v58Transform(b,point)))<1e-10,'skin matrix matches CPU transform');
 const fixture=new ArrayBuffer(134),dv=new DataView(fixture);dv.setUint32(80,1,true);
 [[-14,-72,1341],[-24,-72,1341],[-14,-82,1341]].forEach((v,i)=>v.forEach((n,j)=>dv.setFloat32(96+i*12+j*4,n,true)));
 const parsed=v58Parse(fixture,{...spec,bytes:134});
 ok(parsed.data.faces.length===3&&parsed.weights.every(x=>x>=0&&x<=1),'STL cluster and bounded weights');
 await v58Load(async e=>{if(e===V58_DEEP[3])throw Error('simulated missing');return parsed});
 ok(v58Meshes.size===16&&v58Failed.length===1,'partial deep load independent of skeleton');
 let retry=0;await v58Load(async()=>{retry++;return parsed});
 ok(retry===1&&v58Meshes.size===17,'retry only missing muscle');
 $('v58Subscap').checked=true;settings.cuffOnly=true;ok(V58_DEEP.filter(v58Visible).length===1,'cuff isolation includes subscapularis');settings.cuffOnly=false;
 v572Preset('neutral');state.yaw=Math.PI/2;
 const drag={accum:0,planeAccum:0};
 ok(v58DragShoulder(drag,100,-100)&&Number($('elev').value)===43&&Number($('plane').value)===43,'diagonal shoulder drag changes both axes');
 v58DragShoulder(drag,1000,-1000);ok(Number($('plane').value)===90&&Number($('elev').value)===180,'shoulder drag clamps both axes');
 v58DragShoulder(drag,-1000,1000);ok(Number($('plane').value)===0&&Number($('elev').value)===0,'reverse drag returns to neutral');
 setDisplayPreset('skeleton');ok(V58_DEEP.every(s=>!$(s.id).checked),'skeleton hides all deep muscles');
 }

 {
 ok(V59_ALL.length===9&&new Set(V59_ALL.map(s=>s.file)).size===9,'nine distinct atlas muscle regions');
 const idx=v59PointIndex([[1,2,3],[4,5,6]]);
 ok(idx.near([1,2,3])===0&&idx.near([20,20,20])===Infinity,'spatial attachment index');
 const s=V58_DEEP[3],box=[[0,0,0],[2,2,2]],p=[1,1,1];
 const c={ribs:v59PointIndex([p]),surface:null,origin:null,insertion:[[2,1,1],[2,2,1]]};
 ok(v59Weight(s,p,box,c)===0,'rib-near origin held stationary');
 ok(v59Weight(s,[2,1,1],box,c)===1,'scapular candidate held to moving anchor');
 const rest=[1,1,0],posed=[1.5,1,0],center=[0,0];
 ok(len(sub(v59Wrap(rest,posed,0,center),posed))<1e-9&&len(sub(v59Wrap(rest,posed,1,center),posed))<1e-9,'radial guard preserves exact endpoints');
 ok(len(sub(v59Wrap(rest,rest,.5,center),rest))<1e-9,'radial guard preserves neutral');
 ok(Math.abs(v59Wrap(rest,posed,.5,center)[0]-1.067)<1e-9,'mid-belly outward radial movement capped');
 ok(Math.abs(v59Wrap(rest,[.5,1,0],.5,center)[0]-.972)<1e-9,'mid-belly inward radial movement capped');
 ok(V59_SURFACE.every(s=>v58Meshes.has(s.id)),'surface meshes independently loaded');
 $('v56Back').checked=true;$('v59TrapUpper').checked=true;
 ok(v58Visible(V59_SURFACE[0]),'surface master and individual visibility');
 $('v56Back').checked=false;ok(!v58Visible(V59_SURFACE[0]),'surface master hides atlas muscle');
 setDisplayPreset('skeleton');ok(V59_ALL.every(s=>!$(s.id).checked),'skeleton hides all nine muscle regions');
 }

 {
 ok(V60_ALL.length===17&&new Set(V60_ALL.map(s=>s.file)).size===17,'17 unique atlas regions');
 ok(V60_CHEST.length===8&&V60_DELT.length===3,'eight chest/shoulder regions');
 ok(V60_CHEST.every(s=>v58Meshes.has(s.id)),'all new synthetic meshes loaded');
 const spec=V60_CHEST[0],c=v60WeightContext(spec);
 ok(v60Weight(spec,c.origin[0],null,c)===0&&v60Weight(spec,c.insertion[0],null,c)===1,'chest origin/insertion weights');
 for(const s of V60_CHEST){const c=v60WeightContext(s),w=v60Weight(s,[0,1,0],null,c);ok(Number.isFinite(w)&&w>=0&&w<=1,'finite weights '+s.id)}
 const a={q:ID,p:OR,rest:OR},b={q:qa([0,0,1],1.1),p:[1,2,3],rest:OR},v=[[0,0,0],[1,2,3],[2,1,0],[.5,1,1]],w=[0,.23,.65,1];
 const pts=v60BoundPoints(v60WeightBounds(v,w),a,b),box=v58Box(pts);
 ok(v.every((p,i)=>v58Blend(a,b,p,w[i]).every((x,j)=>x>=box[0][j]-1e-9&&x<=box[1][j]+1e-9)),'weight interval boxes enclose skinned vertices');
 const before=$('muscles').checked;v60Layer('surface');
 ok($('v60PecClav').checked&&$('v60DeltAnt').checked&&!$('v60PecMinor').checked&&!$('v58Serratus').checked,'surface layer excludes deep');
 v60Layer('deep');ok(!$('v60PecClav').checked&&!$('deltoid').checked&&$('v60PecMinor').checked&&$('v58Subscap').checked,'deep layer excludes surface');
 v60Layer('all');ok(V60_ALL.every(s=>$(s.id).checked)&&$('muscles').checked===before,'all layers preserve upper arm switch');
 setDisplayPreset('skeleton');ok(V60_ALL.every(s=>!$(s.id).checked),'skeleton hides all atlas muscles');
 ok(!v60ScapAudit().ready,'missing thorax not claimed as measured clearance');
 ok(!v60DeltReady(),'software harness cannot claim atlas deltoid GPU rendering');
 }

 {
 const v=[[0,0,0],[1,0,0],[0,1,0],[1,1,0]],f=[0,1,2,1,2,3],raw=[0,.15,.85,1];
 const w=v601SmoothWeights(v,f,raw);
 ok(w[0]===0&&w[3]===1,'weight smoothing preserves attachment endpoints');
 ok(w[1]>.15&&w[2]<.85&&Array.from(w).every(x=>x>=0&&x<=1),'smoothing reduces local jumps and stays bounded');
 ok(Array.from(v601SmoothWeights(v,f,[.4,.4,.4,.4])).every(x=>Math.abs(x-.4)<1e-6),'constant field preserved');
 const spec=V58_DEEP[3],p=[1,1,1],c={ribs:v59PointIndex([p]),surface:null,origin:null,insertion:[[2,1,1],[2,2,1]]};
 ok(v601Weight(spec,p,null,c,.5)===0,'serratus nearby rib fixed');
 ok(v601Weight(spec,[2,1,1],null,c,.5)===1,'serratus insertion fully follows scapula');
 const q=[1+atlasScale*16,1,1],withSample=v601Weight(spec,q,null,c,.5),noSample=v601Weight(spec,q,null,{...c,ribs:v59PointIndex([])},.5);
 ok(Math.abs(withSample-noSample)<1e-9,'rib weight finite support continuous at search cutoff');
 const scap=V58_DEEP[4],c2={...c,ribs:v59PointIndex([]),surface:v59PointIndex([p])};
 ok(v601Weight(scap,p,null,c2,.5)===0,'subscapularis near scapula remains attached');
 const pec=V60_CHEST[0],pc=v60WeightContext(pec);
 ok(v601Weight(pec,p,null,pc,0)===0&&v601Weight(pec,p,null,pc,1)===1&&v601Weight(pec,p,null,pc,.5)<.5,'pectoral origin bias preserves endpoints');
 const samples=v601EdgeSamples(v,f,raw,w,10);
 ok(samples.length===5&&samples.every(e=>e.rest>0),'edge audit deduplicates sample edges');
 const a={q:ID,p:OR,rest:OR};
 ok(samples.every(e=>Math.abs(v601EdgeRatio(e,a,a,e.w)-1)<1e-9),'neutral deformation preserves edge lengths');
 const b={q:qa([0,1,0],1.1),p:[1,2,3],rest:OR};
 ok(samples.every(e=>Number.isFinite(v601EdgeRatio(e,a,b,e.w,[0,0]))),'finite guarded edge ratios');
 const s=v601Stats([.1,.9,1,1.2,3]);ok(s.n===5&&s.outside===2&&s.p95===3,'edge statistics distinguish heuristic outliers');
 ok(v601Stats([])===null,'missing edge samples not treated as zero');
 const timing=v601PerfStats([1,2,3,4,10]);ok(timing.n===5&&timing.median===3&&timing.p95===10,'CPU render timing summary');
 ok(v601PerfStats([])===null,'empty timing not measured');
 }

 {
 const sphere=[];for(let i=0;i<20;i++)for(let j=0;j<12;j++){const a=i*Math.PI/10,b=j*Math.PI/6;sphere.push([Math.cos(a)*Math.sin(b),2*Math.cos(b),Math.sin(a)*Math.sin(b)])}
 const ell=v61FitThorax(sphere);
 ok(ell&&Math.abs(v61Radial([1,0,0],ell))<1e-8,'thorax proxy surface radial residual');
 ok(v61Radial([.5,0,0],ell)<0&&v61Radial([1.5,0,0],ell)>0,'signed inward/outward proxy');
 ok(v61FitThorax([])===null,'missing samples never produce a guide');
 const ac=[.8,0,0],pts=[[0,0,.4],[0,.2,.4],[-.2,.3,.5],[-.3,-.2,.6],[-.2,.1,.5],[-.15,-.4,.55],[-.1,.4,.4],[-.3,.4,.3]];
 const data={e:ell,v:pts.map(p=>({p,rest:v61Radial(add(ac,p),ell)}))};
 const q=qa([0,1,0],.5),g=v61Guide(q,ac,120,data);
 ok(g.after<=g.before+1e-9&&g.angles.every(a=>Math.abs(a)<=8),'guide objective non-increasing and correction bounded');
 ok(len(sub(transform(g.q,ac,OR),ac))<1e-9,'guide preserves AC pivot');
 ok(v61Guide(ID,ac,0,data).q===ID,'guide preserves neutral pose');
 ok(v61Guide(q,ac,120,null).q===q,'missing geometry retains old kinematics');
 const v=[[0,0,0],[.1,0,0],[.2,0,0],[.1,.1,0]],f=[0,1,3,1,2,3],old=[0,.1,1,.8],refined=[0,.4,1,.6];
 const s=v61SelectWeights(V60_CHEST[0],v,f,old,refined);
 ok(s.report&&s.report.score<=s.report.baseline+1e-9,'sampled multi-pose choice never worse than baseline objective');
 ok(s.weights[0]===0&&s.weights[2]===1&&Array.from(s.weights).every(w=>w>=0&&w<=1),'selection preserves exact attachment candidates');
 ok(v61Objective([1,1])===0&&v61Objective([.5,2])>0,'symmetric strain objective');
 ok(Math.abs(v61Objective([.5])-v61Objective([2]))<1e-9,'compression and extension scored symmetrically');
 ok(v61EdgeIndices(v,f).length===5,'unique sample edge selection');
 const report=v61AuditReport();ok(report.version==='v0.62'&&!report.collisionValidated&&!report.volumeValidated,'audit declares validation limits');
 }

 {
 const a={q:ID,p:OR,rest:OR},b={q:qa([0,1,0],Math.PI),p:[.1,.2,.3],rest:[.3,.2,.1]},p=[.4,.5,.6],q=[.8,.5,.6],close=(a,b)=>len(sub(a,b))<1e-8;
 ok(close(v62DQBlend(a,b,p,0),v58Transform(a,p))&&close(v62DQBlend(a,b,p,1),v58Transform(b,p)),'DQS attachment endpoints');
 ok(close(v62DQBlend(a,a,p,.43),p),'DQS neutral preserved');
 const neg={...b,q:b.q.map(x=>-x)};
 ok(close(v62DQBlend(a,b,p,.5),v62DQBlend(a,neg,p,.5)),'DQS antipodal invariance');
 ok(Math.abs(len(sub(v62DQBlend(a,b,p,.5),v62DQBlend(a,b,q,.5)))-len(sub(p,q)))<1e-8,'equal-weight DQS preserves local length under twist');
 ok(len(sub(v58Blend(a,b,p,.5),v58Blend(a,b,q,.5)))<1e-8,'linear blend twist collapse fixture');
 const route={kind:1,rest:OR,radii:[1,1]},outside=[2,0,0];
 ok(close(v62RoutePoint([1,0,0],outside,0,route),outside)&&close(v62RoutePoint([1,0,0],outside,1,route),outside),'rib route preserves fixed endpoints');
 const fixed=v62RoutePoint([1,0,0],outside,.5,route);
 ok(fixed[0]<outside[0]&&len(sub(fixed,outside))<=.120000001,'rib route bounded outward correction');
 ok(close(v62RoutePoint(p,p,.5,route),p),'rib route neutral unchanged');
 const r={kind:2,rest:OR,now:OR,radius:1},inside=[.2,0,0],result=v62RoutePoint([.8,0,0],inside,.5,r);
 ok(result[0]>inside[0]&&len(sub(result,inside))<=.120000001,'head proxy bounded inward correction');
 ok(close(v62RoutePoint(p,p,.5,r),p),'head proxy neutral unchanged');
 const vertices=[p,q,[.1,.2,.3],[.1,-.3,.5]],w=[.13,.22,.48,.92],bins=v60WeightBounds(vertices,w),bb=v58Box(v62Bounds(bins,a,b,.12));
 ok(vertices.every((p,i)=>v62DQBlend(a,b,p,w[i]).every((x,j)=>x>=bb[0][j]-1e-9&&x<=bb[1][j]+1e-9)),'DQS camera bounds enclose transformed points');
 const edges=v601EdgeSamples(vertices,[0,1,2,1,2,3],w,w,12),s=v62SelectMode(V60_CHEST[0],edges);
 ok(s&&s.on.score<=s.on.baseline+1e-9&&s.off.score<=s.off.baseline+1e-9,'auto nonlinear candidates include linear baseline');
 $('v62Mode').value='linear';ok(!v62Mode(V60_CHEST[0],{nonlinear:s}).dual&&!v62Mode(V60_CHEST[0],{nonlinear:s}).route,'linear rollback mode');
 $('v62Mode').value='dual';ok(v62Mode(V60_CHEST[0],{nonlinear:s}).dual,'manual DQS mode');
 $('v62Mode').value='auto';ok(v62Mode(V60_CHEST[0],{nonlinear:s}).dual===s[$('v62Route').checked?'on':'off'].dual,'auto mode matches selected result');
 ok(!v62Audit().collisionValidated&&!v62Audit().volumeValidated,'nonlinear audit declares limitations');
 }
 return passed;
};`;
js=js.replace(/\}\)\(\);\s*$/,instrument+"\n})();");
new Function("document","window","fetch","requestAnimationFrame","matchMedia","innerWidth","innerHeight","setTimeout","clearTimeout","console",js)(doc,win,async()=>{throw Error("offline test")},f=>callbacks.push(f),()=>({matches:true}),390,740,()=>0,()=>{}, {log(){},warn(){},error(...a){errors.push(a.map(String).join(" "))}});
assert.equal(errors.length,0,errors.join("\n"));
const results=await win.__runV56Tests();
console.log('PASS v0.62: '+results.length+' integration checks — '+results.join('; '));
