// Run from repository root: node tests/scapular-muscles-v051.mjs
// Checks mesh continuity and moving-bone attachment under SYNTHETIC anatomy.
// This is NOT a validation against real human attachments or bone collisions.
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
const html=readFileSync('index.html','utf8');
const js=html.match(/<script>([\s\S]*?)<\/script>/)?.[1];
assert.ok(js,'inline app script not found');new Function(js);
const groupSource=js.slice(js.indexOf('const SCAPULAR_SPEC='),js.indexOf('const MUSCLE_PARTS='));
const spec=new Function(groupSource+'return {SCAPULAR_SPEC,SCAPULAR_CUFF,SCAPULAR_MUSCLES,SCAPULAR_TENDONS,SCAPULAR_PARTS}')();
assert.equal(spec.SCAPULAR_SPEC.length,4);
assert.deepEqual([...spec.SCAPULAR_CUFF],['supraspinatus','infraspinatus','teresMinor']);
assert.equal(spec.SCAPULAR_SPEC.find(x=>x.id==='teresMajor').cuff,false);
const add=(a,b)=>a.map((x,i)=>x+b[i]),sub=(a,b)=>a.map((x,i)=>x-b[i]);
const mul=(a,k)=>a.map(x=>x*k),dot=(a,b)=>a.reduce((v,x,i)=>v+x*b[i],0);
const len=a=>Math.hypot(...a),norm=a=>mul(a,1/Math.max(len(a),1e-12));
const mix=(a,b,t)=>a+(b-a)*t,clamp=(x,a,b)=>Math.min(b,Math.max(a,x));
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const qa=(axis,r)=>[...mul(norm(axis),Math.sin(r/2)),Math.cos(r/2)];
const qr=(q,v)=>{const u=q.slice(0,3),t=mul(cross(u,v),2);return add(v,add(mul(t,q[3]),cross(u,t)))};
const transform=(q,p,v)=>add(p,qr(q,v));
const atlasScale=.00555;
const atlasToWorld=p=>[-(p[0]+14)*atlasScale,1.70+(p[2]-1341)*atlasScale,(p[1]+72)*atlasScale+.08];
const realAnchors={AC:atlasToWorld([-137,-92,1342]),GH:atlasToWorld([-172,-72,1320])};
const ID=[0,0,0,1],realPose={scapula:[ID,realAnchors.AC],humerus:[ID,realAnchors.GH]};
const scalp=[];
for(let x=-210;x<=-20;x+=5)for(let z=1120;z<=1400;z+=5)for(let y of [-95,-70,-48])
 scalp.push(sub(atlasToWorld([x,y,z]),realAnchors.AC));
const realSurfaceVertices={scapula:scalp};const realSurfaceFaces={scapula:[]};
const attachedAtlas=(name,xyz)=>transform(realPose[name][0],realPose[name][1],
 sub(atlasToWorld(xyz),name==='scapula'?realAnchors.AC:realAnchors.GH));
const geomDynamic={};const settings={scapularMuscles:true};const SOFTWARE=true,gl={};
const mesh=(vertices,faces)=>{
 assert.ok(vertices.length>4&&faces.length>0,'no triangles');
 for(const p of vertices)assert.ok(p.length===3&&p.every(Number.isFinite),'non-finite point');
 for(const f of faces)assert.ok(Number.isInteger(f)&&f>=0&&f<vertices.length,'invalid triangle index');
 return {vertices,faces,normals:vertices.map(()=>[0,0,1])};
};
const env={...spec,realSurfaceVertices,realSurfaceFaces,realAnchors,realPose,atlasScale,atlasToWorld,
 add,sub,mul,dot,len,norm,mix,clamp,qr,transform,mesh,settings,SOFTWARE,gl,geomDynamic,
 attachedAtlas,gdynamic:(name,geom)=>geomDynamic[name]=geom,realReady:true};
const start=js.indexOf('const scapularSurfaceCache=');
const end=js.indexOf('function updateRealMuscles(){',start);
assert.ok(start>0&&end>start,'scapular mesh generator missing');
const code=js.slice(start,end);
const api=new Function(...Object.keys(env),code+
 'return {updateScapularMuscles,get:()=>scapularDiagnostics}')(...Object.values(env));
let cases=0;
for(const elev of [0,30,60,90,120,150,180])for(const axial of [-45,0,45]){
 const r=elev*Math.PI/180;
 realPose.scapula=[qa([1,0,0],r/3),realAnchors.AC];
 realPose.humerus=[qa([0,0,1],r+axial*Math.PI/180),realAnchors.GH];
 api.updateScapularMuscles();
 const d=api.get();
 assert.equal(d.muscles.length,4);
 assert.ok(d.maximumJointGap<1e-10,'muscle-to-tendon seam gap');
 for(const name of [...spec.SCAPULAR_MUSCLES,...spec.SCAPULAR_TENDONS])
  assert.ok(geomDynamic[name]?.vertices?.length>0,'missing mesh '+name);
 assert.equal(d.verified,false,'do not claim clinical/anatomical validation');
 cases++;
}
settings.scapularMuscles=false;api.updateScapularMuscles();
for(const name of [...spec.SCAPULAR_MUSCLES,...spec.SCAPULAR_TENDONS])
 assert.equal(geomDynamic[name],undefined,'stale mesh after OFF: '+name);
for(const id of ['scapularMuscles','cuffOnly','scapulaDisplay','scapularBackView',
 'scapularFrontView','scapularMuscleStatus'])
 assert.ok(html.includes('id="'+id+'"'),'missing UI control '+id);
assert.ok(html.includes('上肢3D解剖ラボ · v0.54'),'wrong app version');
assert.ok(js.includes("settings.scapulaDisplay==='translucent'"),'scapular translucency missing');
assert.ok(js.includes('SCAPULAR_CUFF'),'rotator cuff classification missing');
console.log('PASS v0.54: '+cases+' synthetic shoulder poses × 4 muscles; mesh seams; cleanup; cuff distinction; 6 UI controls. No anatomical accuracy certified.');
