import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
for(const p of ['engine/math.js','anatomy/shoulder-rig-data.js','engine/joint-frames.js','engine/pose-field.js','anatomy/shoulder-pose-field-v076.js','engine/shoulder-kinematics.js'])await import('../'+p);
const fixture=JSON.parse(readFileSync(new URL('./fixtures/shoulder-v076-field-calibration.json',import.meta.url))),M=AnatomyMath,K=ShoulderKinematics;
const rig=new K.ShoulderModel(fixture),field=ShoulderPoseFieldData;
assert.equal(createHash('sha256').update(readFileSync(new URL('../engine/shoulder-kinematics.js',import.meta.url))).digest('hex'),field.referenceSourceHash,'Regenerate the pose field after changing its reference solver');
assert.equal(rig.fieldReason,'validated-field');
assert.equal(rig.fieldSignature,field.signature);
assert.equal(field.values.length,61*19*6);
assert(field.validation.samples===9720&&field.validation.maxRigidErrorBoundMm<1);
const error=(a,b)=>M.len(M.sub(a,b))/rig.scale;
let max=0;
for(const [elev,plane] of [[0,0],[6.75,2.5],[42.75,41.25],[90,40],[120.75,86.25],[130.5,88.75],[150,90],[179.25,88.75],[180,90]]){
 const fast=rig.pose({elev,plane}),ref=rig.pose({elev,plane,backend:'reference'});
 assert.equal(fast.diagnostics.backend,'validated-field');
 assert.equal(ref.diagnostics.backend,'reference');
 for(const p of [fast,ref]){
  assert(error(M.apply(p.bones.clavicle,rig.landmarks.AC),M.apply(p.bones.scapula,rig.landmarks.AC))<1e-8);
  assert(error(M.apply(p.bones.scapula,rig.landmarks.GH),M.apply(p.bones.humerus,rig.landmarks.GH))<1e-8);
 }
 max=Math.max(max,error(fast.AC,ref.AC),error(fast.GH,ref.GH));
 assert(error(fast.AC,ref.AC)<1&&error(fast.GH,ref.GH)<1);
}
// A changed source calibration must use the reference solve, never silently
// apply a pose field calibrated for different bones/contact samples.
const changed=structuredClone(fixture);changed.landmarks.SC[0]+=.001;
const fallback=new K.ShoulderModel(changed);
assert.equal(fallback.fieldReason,'calibration-signature-mismatch');
assert.equal(fallback.pose({elev:90,plane:40}).diagnostics.backend,'reference');
assert.equal(rig.pose({elev:90,plane:40,fit:false}).diagnostics.backend,'reference');
// C1 interpolation: derivatives from both sides of a cell boundary agree.
for(const [e,p,axis]of [[60,42.5,0],[61.5,40,1],[120,87.5,0]]){
 const h=.0001,at=d=>ShoulderPoseField.sample(field,e+(axis===0?d:0),p+(axis===1?d:0));
 const a=at(-h),b=at(0),c=at(h);
 for(let j=0;j<6;j++)assert(Math.abs((b[j]-a[j])/h-(c[j]-b[j])/h)<.001,'joint derivative discontinuity');
}
console.log(JSON.stringify({fieldNodes:61*19,validationSamples:9720,maxRigidErrorBoundMm:field.validation.maxRigidErrorBoundMm,maxSampleJointErrorMm:max,clinicalValidation:false}));
