import assert from 'node:assert/strict';
import {readFileSync}from'node:fs';
for(const f of ['engine/math','anatomy/shoulder-rig-data','engine/joint-frames','anatomy/thorax-contact-data-v076','engine/convex-contact','engine/shoulder-kinematics'])await import('../'+f+'.js');
const fixture=JSON.parse(readFileSync(new URL('./fixtures/shoulder-v076-field-calibration.json',import.meta.url))),rig=new ShoulderKinematics.ShoulderModel(fixture),C=ConvexContact,M=AnatomyMath;
assert(rig.body&&rig.bodySamples.length>900);
const gaps=rig.body.points.map(p=>C.gapMm(p,rig.body,rig.scale));
assert(Math.min(...gaps)>0,'Source neutral placement must already be feasible');
const body=rig.body,dup={...body,normals:[...body.normals,...body.normals],supportsMm:[...body.supportsMm,...body.supportsMm]};
let maxGradientError=0;
for(let i=0;i<150;i++){
 const p=[-.5+(i%13)*.13,.6+(i%17)*.09,-.3+(i%19)*.07],full=C.evaluate(p,body,rig.scale,true),fast=C.penetration(p,body,rig.scale,true);
 assert(Math.abs(Math.min(0,full.gap)-Math.min(0,fast.gap))<1e-10,'Broad phase must preserve the exact unilateral residual');
 assert(Math.abs(C.gapMm(p,dup,rig.scale)-full.gap)<1e-10,'Uniform halfspace duplication invariant');
 assert(M.len(full.gradient)<=1+1e-12);
 for(let j=0;j<3;j++){const lo=p.slice(),hi=p.slice(),h=.001;lo[j]-=h*rig.scale;hi[j]+=h*rig.scale;const d=(C.gapMm(hi,body,rig.scale)-C.gapMm(lo,body,rig.scale))/(2*h);maxGradientError=Math.max(maxGradientError,Math.abs(d-full.gradient[j]));}
}
assert(maxGradientError<1e-5,'Analytic spatial derivative versus centred mm finite differences');
const params=ShoulderKinematics.prior(120,40),{r,J}=rig.jacobian(params,params);
let maxJacobianError=0;
for(let j=0;j<6;j++){const a=params.slice(),b=params.slice(),h=.0001;a[j]-=h;b[j]+=h;const lo=rig.residual(a,params),hi=rig.residual(b,params);for(let i=0;i<r.length;i++)maxJacobianError=Math.max(maxJacobianError,Math.abs((hi[i]-lo[i])/(2*h)-J[j][i]));}
assert(maxJacobianError<1e-4,'Rigid-transform derivative and support gradient must reproduce the full residual');
// Captured interpolation failure interval: converged model must not retain the
// former 24-iteration energy-increasing update branch.
const fine=new ShoulderKinematics.ShoulderModel({...fixture,solverIterations:96});
let maxPositionDifferenceMm=0;
for(const elev of [135,135.75,136.5]){const input={elev,plane:8.75},a=rig.pose(input),b=fine.pose(input);maxPositionDifferenceMm=Math.max(maxPositionDifferenceMm,M.len(M.sub(a.AC,b.AC))/rig.scale,M.len(M.sub(a.GH,b.GH))/rig.scale);}
assert(maxPositionDifferenceMm<1,'48 versus 96 iterations; unchanged 1mm engineering tolerance');
console.log(JSON.stringify({neutralMinSupportMm:Math.min(...gaps),samples:gaps.length,maxGradientError,maxJacobianError,maxPositionDifferenceMm,actualSurfaceCollision:'separate source-triangle audit; support is not Euclidean gap'}));
