// Numerical convergence/performance audit, not a phone FPS or clinical test.
import { readFileSync,writeFileSync } from 'node:fs';
for(const p of ['engine/math.js','anatomy/shoulder-rig-data.js','engine/joint-frames.js','engine/pose-field.js','anatomy/shoulder-pose-field-v076.js','anatomy/thorax-contact-data-v076.js','engine/convex-contact.js','engine/shoulder-kinematics.js'])await import('../../'+p);
const fixture=JSON.parse(readFileSync(process.argv[2]||'tests/fixtures/shoulder-v076-field-calibration.json')),M=AnatomyMath,K=ShoulderKinematics;
const rig=new K.ShoulderModel(fixture),fine=new K.ShoulderModel({...fixture,solverIterations:96}),inputs=[];
for(const plane of [0,40,90])for(let elev=0;elev<=180;elev+=3)inputs.push({elev,plane});
inputs.push({elev:130.5,plane:88.75});
const angle=(a,b)=>2*Math.acos(Math.min(1,Math.abs(a.reduce((s,x,i)=>s+x*b[i],0)))),radius=(lo,hi,p)=>Math.hypot(...p.map((x,i)=>Math.max(Math.abs(x-lo[i]),Math.abs(x-hi[i]))));
const rc=radius([-150.033,-146.766,1331.76],[-11.4439,-47.2239,1366.86],ShoulderRigData.landmarks.SC),rs=radius([-168.685,-99.0763,1201.5],[-59.3178,11.3517,1359.85],ShoulderRigData.landmarks.AC);
let max=0,worst=null;
for(const input of inputs){
 const a=rig.pose({...input,backend:'reference'}),b=fine.pose({...input,backend:'reference'}),ac=M.len(M.sub(a.AC,b.AC))/fixture.scale,gh=M.len(M.sub(a.GH,b.GH))/fixture.scale;
 const bound=Math.max(2*rc*Math.sin(angle(a.qC,b.qC)/2),ac+2*rs*Math.sin(angle(a.qS,b.qS)/2),gh);
 if(bound>max){max=bound;worst={input,acMm:ac,ghMm:gh,boundMm:bound}}
}
const times={};
for(const backend of ['auto','reference']){
 rig.cache.clear();rig.girdleCache.clear();const ms=[];
 for(const input of inputs){const t=performance.now();rig.pose({...input,backend});ms.push(performance.now()-t)}
 ms.sort((a,b)=>a-b);times[backend]={samples:ms.length,medianMs:ms[Math.floor(ms.length/2)],p95Ms:ms[Math.floor(ms.length*.95)],maxMs:ms.at(-1)};
}
const result={samples:inputs.length,iterations:[48,96],maxRigidDifferenceBoundMm:max,worst,times,environment:{node:process.version,physicalPhone:false,renderingMeasured:false},clinicalValidation:false};
writeFileSync('docs/validation/v076/reference-convergence.json',JSON.stringify(result,null,2)+'\n');console.log(result);
if(max>=1)throw Error('48 vs 96 iterations exceeds 1mm numerical tolerance; investigate without relaxing the threshold');
