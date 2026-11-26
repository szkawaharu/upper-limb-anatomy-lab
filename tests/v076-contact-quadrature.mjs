import assert from 'node:assert/strict';
import {readFileSync}from'node:fs';
for(const p of ['engine/math.js','anatomy/shoulder-rig-data.js','engine/joint-frames.js','engine/shoulder-kinematics.js','engine/rig-diagnostics.js'])await import('../'+p);
const fixture=JSON.parse(readFileSync(new URL('./fixtures/shoulder-v076-field-calibration.json',import.meta.url))),K=ShoulderKinematics,M=AnatomyMath,F=JointFrames;
const base=new K.ShoulderModel(fixture),copies=[2,4].map(n=>new K.ShoulderModel({...fixture,contact:{...fixture.contact,points:Array.from({length:n},()=>fixture.contact.points).flat()}}));
const error=(a,b)=>M.len(M.sub(a,b))/fixture.scale;
let maxDensityErrorMm=0,maxHTErrorDeg=0,minClearance=Infinity;
for(const plane of[0,40,90])for(const elev of[0,30,60,90,120,150,180]){
 const input={elev,plane,axial:30},p=base.pose({...input,backend:'reference'}),report=RigDiagnostics.shoulderMotion(base,p,input);
 assert.equal(report.clinicalROMorDiseaseClassification,null);
 assert.equal(report.contributionsAreNotAdditive,true);
 maxHTErrorDeg=Math.max(maxHTErrorDeg,Math.abs(report.humerothoracicElevation-elev));
 assert(Math.abs(report.humerothoracicElevation-elev)<1e-5);
 const st=report.scapulothoracic,angles=[st.internalRotation,-st.upwardRotation,st.posteriorTilt].map(x=>x*Math.PI/180),q=M.qmul(base.thoraxFrame,F.fromYXZ(angles)),actual=M.qmul(p.qS,base.scapulaFrame);
 for(const v of[[1,0,0],[0,1,0],[0,0,1]])assert(M.len(M.sub(M.qrot(q,v),M.qrot(actual,v)))<1e-12,'reported ST angles reconstruct the actual rigid frame');
 const prior=report.unconstrainedPriorScapulothoracic;
 for(const key of['internalRotation','upwardRotation','posteriorTilt'])assert(Math.abs(report.contactChangeFromPrior[key]-(st[key]-prior[key]))<1e-9);
 // Repeating the same spatial quadrature samples supplies no new geometry.
 const r=base.residual(p.parameters,p.prior),loss=r.reduce((s,x)=>s+x*x,0);
 for(const other of copies){
  const rr=other.residual(p.parameters,p.prior),otherLoss=rr.reduce((s,x)=>s+x*x,0);
  assert(Math.abs(loss-otherLoss)<1e-10,'contact objective depends on sample count');
  const x=other.pose({...input,backend:'reference'});
  maxDensityErrorMm=Math.max(maxDensityErrorMm,error(x.AC,p.AC),error(x.GH,p.GH));
  assert(error(x.AC,p.AC)<.0001&&error(x.GH,p.GH)<.0001,'spatial sampling changes the final bones');
 }
 minClearance=Math.min(minClearance,p.diagnostics.absoluteClearanceMm[0]);
}
assert(minClearance>0,'new penetration in the captured proxy samples');
console.log(JSON.stringify({poses:21,copies:[2,4],maxDensityErrorMm,maxHTErrorDeg,minProxyClearanceMm:minClearance,clinicalValidation:false}));
