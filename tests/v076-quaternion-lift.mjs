import assert from 'node:assert/strict';
import {readFileSync}from'node:fs';
for(const file of ['math','fiber-path','volume-constraint','contact-constraint','tendon-deformation','muscle-deformation'])await import('../engine/'+file+'.js');
const M=AnatomyMath,E=MuscleDeformation,P=FiberPath,r=JSON.parse(readFileSync(new URL('./fixtures/v076-quaternion-lift.json',import.meta.url)));
const A=r.from.state,B=r.to.state,x=r.from.sample;
assert(A.a.q.reduce((s,v,i)=>s+v*A.b.q[i],0)>0&&B.a.q.reduce((s,v,i)=>s+v*B.b.q[i],0)<0,'captured hemisphere crossing');
const controls=Array.from({length:A.rest.length/12},(_,i)=>Array.from({length:4},(_,j)=>A.rest.slice(i*12+j*3,i*12+j*3+3)));
const binding={field:{controls},reference:A.reference,sheet:A.sheet,blendRange:A.blendRange,definition:{id:r.muscleId,tissues:{tendonStart:A.tendonStart}}};
const results=[];
for(const divisions of [1,2,4,8,16,32]){
 let previous=null,maxMm=0;
 for(let i=0;i<=divisions;i++){
  const u=i/divisions,b={q:M.qmix(A.b.q,B.b.q,u),p:M.mix(A.b.p,B.b.p,u),rest:A.b.rest},f=E.pose(binding,A.a,b,A.contact),q=E.point(f,x.p,x.t,x.l,x.materialT);
  assert(q.every(Number.isFinite));
  if(previous)maxMm=Math.max(maxMm,M.len(M.sub(q,previous))/.00555);previous=q;
  for(const t of [0,1])assert(M.len(M.sub(E.point(f,x.p,t,x.l,x.materialT),M.apply(t?b:A.a,x.p)))<1e-12,'fixed attachment');
 }
 results.push({divisions,maxMm});
}
assert(results[3].maxMm/results[0].maxMm<.3,'material frame still switches at 180deg');
assert(results.at(-1).maxMm<results[0].maxMm/16,'refinement converges');
// An unspecified full-turn quaternion lift is reported and remains finite;
// this safety path is separate from the ordinary 180deg hemisphere crossing.
const singular=E.pose(binding,M.identity(),{...M.identity(),q:[0,0,0,-1]});
assert(singular.liftSingular);
const q=E.evaluate(singular,x.p,x.t,x.l,null,x.materialT);
assert(q.point.every(Number.isFinite)&&q.invalidFrame);
console.log(JSON.stringify({id:r.muscleId,vertexId:r.vertexId,oldJumpMm:r.stepMm,refinement:results}));
