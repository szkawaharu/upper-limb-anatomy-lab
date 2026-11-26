import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
for (const p of ['engine/math.js','anatomy/shoulder-rig-data.js','engine/joint-frames.js','engine/shoulder-kinematics.js']) await import('../'+p);
const M=globalThis.AnatomyMath,F=globalThis.JointFrames,K=globalThis.ShoulderKinematics,D=globalThis.ShoulderRigData;
const scale=.00555,map=p=>[-(p[0]+14)*scale,1.7+(p[2]-1341)*scale,(p[1]+72)*scale+.08],L=Object.fromEntries(Object.entries(D.landmarks).map(([k,p])=>[k,map(p)]));
const fixture=JSON.parse(readFileSync(new URL('./fixtures/shoulder-v076-atlas.json',import.meta.url))),original=JSON.stringify(fixture);
const gaps=fixture.points.map(p=>K.clearance(p,fixture.surface)/scale).sort((a,b)=>a-b),mid=gaps[Math.floor(gaps.length/2)],dev=gaps.map(x=>Math.abs(x-mid)).sort((a,b)=>a-b),sigma=1.4826*dev[Math.floor(dev.length/2)];
const rig=new K.ShoulderModel({landmarks:L,scale,contact:{...fixture,uncertaintyMm:sigma}});
const err=(a,b)=>M.len(M.sub(a,b))/scale,close=(a,b,t=1e-7)=>assert(err(a,b)<t,JSON.stringify({a,b,mm:err(a,b)}));
// Independent mesh-cap/sphere calibration regressions. These confirm the
// estimated source landmarks, not clinical millimetre accuracy.
assert(M.len(M.sub(D.landmarks.SC,[-14.910239325,-129.664673633,1341.374530473]))<.001);
assert(M.len(M.sub(D.landmarks.SC,[-15,-80,1343]))>49);
assert(M.len(M.sub(D.landmarks.AC,[-147.163892473,-62.159424986,1357.282684326]))<.001);
for(const q of [rig.clavicleFrame,rig.scapulaFrame]){
 const axes=[[1,0,0],[0,1,0],[0,0,1]].map(v=>M.qrot(q,v));
 assert(Math.abs(M.dot(axes[0],M.cross(axes[1],axes[2]))-1)<1e-12,'proper anatomical frame');
 const round=F.fromYXZ(F.yxz(q));
 for(const v of [[1,0,0],[0,1,0],[0,0,1]])close(M.qrot(q,v),M.qrot(round,v));
}
// Mesh origins may differ from the joint pivot on both parent and child.
const parent={q:M.qaxis([1,2,3],.7),p:[1,2,3],rest:[-.7,.2,4]},pivot=[.3,2,-1],origin=[2,.1,3],child=F.attach(parent,pivot,M.qaxis([0,1,0],.2),origin);
close(M.apply(parent,pivot),M.apply(child,pivot));
const shaft=M.norm(M.sub(L.EL,L.GH)),length=M.len(M.sub(L.AC,L.SC));
let count=0,maxContactStep=0,minClearance=Infinity,maxIterationStep=0;
for(const plane of [0,40,90]){
 let previous=null;
 for(let elev=0;elev<=180;elev++){
  const p=rig.pose({elev,plane});
  assert([...p.qC,...p.qS,...p.qH,...p.AC,...p.GH].every(Number.isFinite));
  assert(Math.abs(M.len(M.sub(p.AC,p.SC))-length)<1e-12,'rigid clavicle');
  close(M.apply(p.bones.clavicle,L.SC),L.SC);
  close(M.apply(p.bones.clavicle,L.AC),M.apply(p.bones.scapula,L.AC));
  close(M.apply(p.bones.scapula,L.GH),M.apply(p.bones.humerus,L.GH));
  const e=elev*Math.PI/180,f=plane*Math.PI/180;
  close(M.qrot(p.qH,shaft),M.qrot(rig.thoraxFrame,[Math.sin(e)*Math.sin(f),-Math.cos(e),Math.sin(e)*Math.cos(f)]));
  for(const axial of [-90,90])close(M.qrot(rig.pose({elev,plane,axial,fit:false}).qH,shaft),M.qrot(p.qH,shaft));
  if(previous)maxContactStep=Math.max(maxContactStep,err(p.AC,previous.AC));
  previous=p;
  minClearance=Math.min(minClearance,p.diagnostics.absoluteClearanceMm[0]);
  maxIterationStep=Math.max(maxIterationStep,p.diagnostics.lastIterationStepDeg);
  assert.equal(p.diagnostics.clinicalValidation,false);
  count++;
 }
}
// A pose must not depend on how it was reached, including full inverse return.
const reference=rig.pose({elev:120,plane:40});
rig.cache.clear();for(let elev=180;elev>=0;elev-=3)rig.pose({elev,plane:90});
assert.deepEqual(rig.pose({elev:120,plane:40}),reference);
assert.equal(JSON.stringify(fixture),original,'source geometry remains unchanged');
const zero=rig.pose({elev:0,plane:0});close(zero.AC,L.AC);close(zero.GH,L.GH);
assert(zero.diagnostics.atlasNeutralHumerusOffsetDeg>10,'source atlas offset is calibrated, not assumed zero');
assert.throws(()=>rig.pose({elev:NaN}),/Non-finite/);
// Refining a 3-degree interval must reduce motion. A discrete solution switch
// would retain a finite jump as the sampling interval shrinks.
const refinement=[];
for(const step of [3,1.5,.75,.375]){
 let prior=null,max=0;
 for(let elev=120;elev<=123+1e-9;elev+=step){const p=rig.pose({elev,plane:0});if(prior)max=Math.max(max,err(p.AC,prior.AC));prior=p}
 refinement.push({step,max});
}
assert(refinement.at(-1).max/refinement[0].max<.3,'non-convergent joint motion');
assert(minClearance>0,'new proxy penetration in captured costal sample set');
// Child articulation cannot change any upstream shoulder transform.
for(const elbow of [0,30,90,140])for(const pronation of [-90,0,90]){
 const u=F.attach(reference.bones.humerus,L.EL,M.qaxis(M.norm(M.sub(map(D.elbowAxis[1]),map(D.elbowAxis[0]))),elbow*Math.PI/180));
 const RH=map([-232,-68,1039]),axis=M.norm(M.sub(map([-240,-110,825]),RH)),r=F.attach(u,RH,M.qaxis(axis,pronation*Math.PI/180),L.EL);
 close(M.apply(u,RH),M.apply(r,RH));
 assert.deepEqual(rig.pose({elev:120,plane:40}),reference);
}
console.log(JSON.stringify({poses:count,sourceOffsetDeg:zero.diagnostics.atlasNeutralHumerusOffsetDeg,maxAcStepMmPerDegree:maxContactStep,minProxyClearanceMm:minClearance,maxLastIterationStepDeg:maxIterationStep,refinement,clinicalValidation:false}));
