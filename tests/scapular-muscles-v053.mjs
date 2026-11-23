// Run: node tests/scapular-muscles-v053.mjs
// Rigorous SYNTHETIC test of mesh generation, not anatomical validation.
// Scapular surface = triangulated curved posterior and anterior sheets.
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
const html=readFileSync('index.html','utf8');
const js=html.match(/<script>([\s\S]*?)<\/script>/)?.[1];
assert.ok(js,'inline application script missing');new Function(js);
const group=js.slice(js.indexOf('const SCAPULAR_SPEC='),js.indexOf('const MUSCLE_PARTS='));
const spec=new Function(group+'return{SCAPULAR_SPEC,SCAPULAR_MUSCLES,SCAPULAR_TENDONS,SCAPULAR_CUFF,SCAPULAR_PARTS}')();
assert.deepEqual(spec.SCAPULAR_MUSCLES,['supraspinatus','infraspinatus','teresMinor','teresMajor']);
assert.deepEqual([...spec.SCAPULAR_CUFF],['supraspinatus','infraspinatus','teresMinor']);
const add=(a,b)=>a.map((x,i)=>x+b[i]);
const sub=(a,b)=>a.map((x,i)=>x-b[i]);
const mul=(a,k)=>a.map(x=>x*k);
const dot=(a,b)=>a.reduce((v,x,i)=>v+x*b[i],0);
const len=a=>Math.hypot(...a),norm=a=>mul(a,1/(len(a)||1));
const mix=(a,b,t)=>a+(b-a)*t,clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const qa=(axis,theta)=>[...mul(norm(axis),Math.sin(theta/2)),Math.cos(theta/2)];
const qr=(q,v)=>{const u=q.slice(0,3),t=mul(cross(u,v),2);return add(v,add(mul(t,q[3]),cross(u,t)))};
const qc=q=>[-q[0],-q[1],-q[2],q[3]];
const transform=(q,p,v)=>add(qr(q,v),p);
const atlasScale=.00555;
const atlasToWorld=p=>[-(p[0]+14)*atlasScale,1.7+(p[2]-1341)*atlasScale,(p[1]+72)*atlasScale+.08];
const realAnchors={AC:atlasToWorld([-137,-92,1342]),GH:atlasToWorld([-172,-72,1320])};
const ID=[0,0,0,1],realPose={scapula:[ID,realAnchors.AC],humerus:[ID,realAnchors.GH]};
const realSurfaceVertices={scapula:[]},realSurfaceFaces={scapula:[]};
const posterior=(x,y)=>.15+.020*Math.sin(3*x)*Math.cos(2*y);
const cols=80,rows=90,x0=-.80,dx=.019,y0=-1.36,dy=.020;
for(let layer=0;layer<2;layer++){
 const offset=realSurfaceVertices.scapula.length;
 for(let y=0;y<=rows;y++)for(let x=0;x<=cols;x++){
  const xx=x0+x*dx,yy=y0+y*dy;
  realSurfaceVertices.scapula.push([xx,yy,layer===0?posterior(xx,yy):-.14]);
 }
 for(let y=0;y<rows;y++)for(let x=0;x<cols;x++){
  const a=offset+y*(cols+1)+x,b=a+1,c=a+cols+1,d=c+1;
  realSurfaceFaces.scapula.push(a,b,c,b,d,c);
 }
}
const attachedAtlas=(bone,xyz)=>transform(realPose[bone][0],realPose[bone][1],
 sub(atlasToWorld(xyz),realAnchors[bone==='scapula'?'AC':'GH']));
const geomDynamic={},settings={scapularMuscles:true},SOFTWARE=true,gl={};
const mesh=(vertices,faces)=>{
 assert.ok(vertices.length>5&&faces.length>20,'nonempty generated mesh');
 for(const p of vertices)assert.ok(p.length===3&&p.every(Number.isFinite),'nonfinite point');
 for(const f of faces)assert.ok(Number.isInteger(f)&&f>=0&&f<vertices.length,'invalid triangle index');
 return {vertices,faces};
};
const env={...spec,add,sub,mul,dot,len,norm,mix,clamp,qr,transform,realAnchors,
 atlasScale,atlasToWorld,realSurfaceVertices,realSurfaceFaces,realPose,attachedAtlas,
 mesh,settings,SOFTWARE,gl,geomDynamic,gdynamic:(id,g)=>{geomDynamic[id]=g},realReady:true};
const start=js.indexOf('const scapularSurfaceCache=');
const end=js.indexOf('function updateRealMuscles(){',start);
assert.ok(start>0&&end>start,'v0.53 projector/generator missing');
const segment=js.slice(start,end);
const api=new Function(...Object.keys(env),segment+
 '\nreturn{run:updateScapularMuscles,diagnostics:()=>scapularDiagnostics,cache:scapularSurfaceCache}')(...Object.values(env));
const sampleCount=7,bellyRows=9,surfaceStart=bellyRows*sampleCount,poses=[];
for(const deg of [0,30,60,90,120,150,180])for(const axial of [-45,0,45]){
 const theta=deg*Math.PI/180;
 realPose.scapula=[qa([1,0,0],theta/3),realAnchors.AC];
 realPose.humerus=[qa([0,0,1],theta+axial*Math.PI/180),realAnchors.GH];
 api.run();
 const report=api.diagnostics();
 assert.equal(report.muscles.length,4);
 assert.ok(report.maximumJointGap<1e-9,'muscle-to-tendon seam gap');
 assert.equal(report.verified,false);
 assert.equal(report.projectionFallbacks,0,'test target must hit posterior triangles');
 assert.ok(report.rayProjectionHits>50,'ray projection was not exercised');
 assert.ok(report.minimumPosteriorProjectionClearanceMm>0,'posterior offset nonpositive');
 assert.ok(Number.isInteger(report.approximateHumeralHeadGuideMoves),'missing humeral-head visual guide count');
 let minClearance=Infinity;
 for(const id of spec.SCAPULAR_MUSCLES){
  const g=geomDynamic[id],t=geomDynamic['tendon'+id[0].toUpperCase()+id.slice(1)];
  assert.ok(g&&t,'missing mesh '+id);
  const centerRow=4*sampleCount,backOffset=bellyRows*sampleCount;
  const edgeThickness=len(sub(g.vertices[centerRow],g.vertices[backOffset+centerRow]))/2;
  const middleThickness=len(sub(g.vertices[centerRow+3],g.vertices[backOffset+centerRow+3]))/2;
  assert.ok(middleThickness>edgeThickness*2,'hard rectangular flap edge persists: '+id);
  // Station t=1 must follow the pose of the humerus precisely, despite
  // the approximate spherical head clearance guide for intermediate stations.
  const row=report.muscles.find(m=>m.id===id);
  assert.ok(row?.endTrack?.every(p=>p.every(Number.isFinite)));
  // The bottom half of the belly at scapula-attached rows 0..6 MUST be
  // outside the curved posterior envelope, never behind the bone.
  for(let k=0;k<=6;k++)for(let j=0;j<sampleCount;j++){
   const point=g.vertices[surfaceStart+k*sampleCount+j];
   const local=qr(qc(realPose.scapula[0]),sub(point,realPose.scapula[1]));
   minClearance=Math.min(minClearance,local[2]-posterior(local[0],local[1]));
  }
 }
 assert.ok(minClearance>.001,'red patches: muscle under posterior bone surface');
 poses.push({deg,axial,minClearance});
}
assert.equal(poses.length,21);
settings.scapularMuscles=false;api.run();
for(const id of [...spec.SCAPULAR_MUSCLES,...spec.SCAPULAR_TENDONS])
 assert.equal(geomDynamic[id],undefined,'stale mesh after OFF: '+id);
for(const id of ['scapularMuscles','cuffOnly','scapulaDisplay','scapularBackView','scapularFrontView','scapularMuscleStatus'])
 assert.ok(html.includes('id="'+id+'"'),'missing control '+id);
assert.ok(html.includes('上肢3D解剖ラボ · v0.53'),'v0.53 not displayed');
assert.ok(js.includes("settings.scapulaDisplay==='translucent'"),'scapula opacity mode removed');
console.log('PASS v0.53: 21 synthetic shoulder poses × 4 muscles, positive posterior clearance, 0 projected misses, 0 tendon seams, OFF cleanup and UI.');
