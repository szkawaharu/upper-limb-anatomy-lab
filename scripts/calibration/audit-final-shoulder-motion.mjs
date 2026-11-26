// Independent reporting of final ST angles and their deviation from the prior.
// A zero collision count must not stand in for a plausible anatomical pose.
import{readFileSync,writeFileSync,mkdirSync}from'node:fs';import{createHash}from'node:crypto';
const files=['engine/math','anatomy/shoulder-rig-data','engine/joint-frames','engine/pose-field','anatomy/shoulder-pose-field-v076','anatomy/thorax-contact-data-v076','engine/convex-contact','engine/shoulder-kinematics','engine/rig-diagnostics'];
for(const f of files)await import('../../'+f+'.js');
const input='tests/fixtures/shoulder-v076-field-calibration.json',fixture=JSON.parse(readFileSync(input)),model=new ShoulderKinematics.ShoulderModel(fixture),loftOnly=new ShoulderKinematics.ShoulderModel(fixture),rows=[];
// Offline ablation only: disable the added safeguard to distinguish its effect
// from the original rest-relative loft term. Never used by the application.
loftOnly.body=null;loftOnly.field=null;
for(const plane of[0,40,90])for(const elev of[0,30,60,90,120,150,180]){
 const input={elev,plane},pose=model.pose(input),motion=RigDiagnostics.shoulderMotion(model,pose,input);
 const before=loftOnly.pose({...input,backend:'reference'}),loftMotion=RigDiagnostics.shoulderMotion(loftOnly,before,input);
 rows.push({input,parameters:pose.parameters,prior:pose.prior,parameterDelta:pose.parameters.map((x,j)=>x-pose.prior[j]),motion,loftOnlyAblation:{parameters:before.parameters,motion:loftMotion,sourceBoneCollisionCertified:false},connections:RigDiagnostics.connections(model,pose.bones),backend:pose.diagnostics.backend});
}
const out=process.argv[2]||'test-results/v076-source-bone-motion.json';mkdirSync('test-results',{recursive:true});
const sha=p=>createHash('sha256').update(readFileSync(p)).digest('hex'),result={method:'Final SC/AC/ST/GH numerical coordinate report; source collision safeguards and prior deviation reported separately; no clinical certification',sourceSha256:Object.fromEntries([input,...files.map(x=>x+'.js')].map(p=>[p,sha(p)])),rows};writeFileSync(out,JSON.stringify(result,null,2));
console.log(JSON.stringify({poses:rows.length,maxAbsoluteParameterDeltaDeg:Array.from({length:6},(_,j)=>Math.max(...rows.map(r=>Math.abs(r.parameterDelta[j])))),at180:rows.filter(x=>x.input.elev===180).map(x=>({plane:x.input.plane,ST:x.motion.scapulothoracic,STChangeFromPrior:x.motion.contactChangeFromPrior})),maxSharedPivotErrorMm:Math.max(...rows.flatMap(r=>r.connections.map(c=>c.sharedPivotErrorMm)))}));
