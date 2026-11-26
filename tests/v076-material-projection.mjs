// Diagnostic experiment only; this projection is NOT used by the production binding.
import assert from 'node:assert/strict';
import { closestMaterial } from '../scripts/experimental/material-projection.mjs';
import { readFileSync } from 'node:fs';
for(const p of ['math','fiber-path','volume-constraint','contact-constraint','tendon-deformation','muscle-deformation'])await import('../engine/'+p+'.js');
const d=JSON.parse(readFileSync(new URL('./fixtures/v076-material-projection.json',import.meta.url))),M=AnatomyMath,P=FiberPath,E=MuscleDeformation,x=d.from.sample,s=d.from.state;
const controls=Array.from({length:s.rest.length/12},(_,i)=>Array.from({length:4},(_,j)=>s.rest.slice(i*12+j*3,i*12+j*3+3))),c=P.laneCurve(controls,x.l),t=closestMaterial(c,x.p,x.t),R=P.frame(controls,x.l,t,s.reference);
const residual=M.dot(M.sub(x.p,R.c),R.t)/.00555;
assert(Math.abs(residual)<.001,'closest-point orthogonality lost');
assert(t>.96&&t<.99);
const B={field:{controls},reference:s.reference,sheet:s.sheet,blendRange:s.blendRange,definition:{tissues:{tendonStart:s.tendonStart}}};
const point=(f,material)=>E.point(E.pose(B,f.a,f.b,f.contact),x.p,x.t,x.l,material);
const oldMm=M.len(M.sub(point(d.from.state,x.materialT),point(d.to.state,x.materialT)))/.00555,newMm=M.len(M.sub(point(d.from.state,t),point(d.to.state,t)))/.00555;
// The geometry has not changed or been thinned. Only the rest chart location
// is corrected. Bone weights/Dirichlet ends are exactly the same.
assert(newMm<oldMm/2,'captured material chart amplification remains');
for(const w of [0,1]){
 const a=E.point(E.pose(B,s.a,s.b,s.contact),x.p,w,x.l,t);
 assert(M.len(M.sub(a,M.apply(w?s.b:s.a,x.p)))<1e-12);
}
const line=[[0,0,0],[1/3,0,0],[2/3,0,0],[1,0,0]];
for(const u of [0,.00001,.01,.4,.99999,1])assert(Math.abs(closestMaterial(line,[u,1,0])-u)<1e-6);
const identity=E.pose(B,M.identity(),M.identity());
assert(M.len(M.sub(E.point(identity,x.p,x.t,x.l,t),x.p))<1e-12,'rest anatomy preserved');
console.log(JSON.stringify({id:d.muscleId,vertexId:d.vertexId,oldMaterialT:x.materialT,orthogonalMaterialT:t,tangentialResidualMm:residual,oldStepMm:oldMm,newStepMm:newMm,anatomicalChartValidation:false}));
