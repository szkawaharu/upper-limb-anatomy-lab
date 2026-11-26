// Compare identical geometry under the original sum and corrected spatial mean.
// The reported rotation is kinematics, not a frozen-shoulder diagnosis.
import{readFileSync,writeFileSync,mkdirSync}from'node:fs';
import{execFileSync}from'node:child_process';
import vm from'node:vm';
for(const p of['engine/math.js','anatomy/shoulder-rig-data.js','engine/joint-frames.js','engine/pose-field.js','engine/shoulder-kinematics.js','engine/rig-diagnostics.js'])await import('../../'+p);
const baseRef='6c5ed7ce68aaea64da458e6ba5b58bbbb748df68',fixture=JSON.parse(readFileSync('tests/fixtures/shoulder-v076-field-calibration.json')),scope={AnatomyMath,JointFrames,ShoulderRigData,ShoulderPoseField};
vm.runInNewContext(execFileSync('git',['show',baseRef+':engine/shoulder-kinematics.js'],{encoding:'utf8'}),scope);
const M=AnatomyMath,scale=fixture.scale,rows=[],density=[];
for(const[label,K]of[['before',scope.ShoulderKinematics],['after',ShoulderKinematics]]){
 const model=new K.ShoulderModel(fixture);
 for(const plane of[0,40,90])for(const elev of[0,30,60,90,120,150,180]){
  const input={elev,plane},p=model.pose({...input,backend:'reference'});
  rows.push({model:label,input,motion:RigDiagnostics.shoulderMotion(model,p,input),contact:p.diagnostics,AC:p.AC,GH:p.GH});
 }
 for(const copies of[2,4]){
  const duplicate=new K.ShoulderModel({...fixture,contact:{...fixture.contact,points:Array.from({length:copies},()=>fixture.contact.points).flat()}});
  for(const plane of[0,40,90]){
   const input={elev:120,plane},a=model.pose({...input,backend:'reference'}),b=duplicate.pose({...input,backend:'reference'}),r=RigDiagnostics.shoulderMotion(model,a,input),s=RigDiagnostics.shoulderMotion(duplicate,b,input);
   density.push({model:label,copies,input,ACdifferenceMm:M.len(M.sub(a.AC,b.AC))/scale,GHdifferenceMm:M.len(M.sub(a.GH,b.GH))/scale,STupwardDifferenceDeg:s.scapulothoracic.upwardRotation-r.scapulothoracic.upwardRotation});
  }
 }
}
const out='docs/validation/v076/contact-rms';mkdirSync(out,{recursive:true});
const sweepRig=new ShoulderKinematics.ShoulderModel(fixture),sweep=[];
for(const plane of[0,40,90])for(let elev=0;elev<=180;elev++){
 const p=sweepRig.pose({elev,plane,backend:'reference'}),d=p.diagnostics;
 sweep.push({plane,elev,absoluteClearanceMm:d.absoluteClearanceMm,clearanceDeltaMm:d.clearanceDeltaMm,clearanceRmsMm:d.clearanceRmsMm,extrapolated:d.surfaceExtrapolatedPoints,lastStepDeg:d.lastIterationStepDeg});
}
const proxySweep={poses:sweep.length,minClearanceMm:Math.min(...sweep.map(r=>r.absoluteClearanceMm[0])),maxClearanceMm:Math.max(...sweep.map(r=>r.absoluteClearanceMm[1])),maxOutside:Math.max(...sweep.map(r=>r.extrapolated)),maxLastStepDeg:Math.max(...sweep.map(r=>r.lastStepDeg)),rows:sweep,exactBoneCollisionEvaluated:false};
const result={beforeSource:baseRef,sourceInput:'tests/fixtures/shoulder-v076-field-calibration.json',sampleCount:fixture.contact.points.length,rows,density,proxySweep,photoVersion:'0.75.0',exactPhotoAnglesRecovered:false,clinicalDiseaseClassification:null,forceOrActivationEstimated:false};
writeFileSync(out+'/shoulder-motion.json',JSON.stringify(result,null,2)+'\n');
console.log(JSON.stringify({density,before120:rows.find(x=>x.model==='before'&&x.input.elev===120&&x.input.plane===90).motion,after120:rows.find(x=>x.model==='after'&&x.input.elev===120&&x.input.plane===90).motion},null,2));
