// Offline generation. Run after an actual-atlas browser export, then commit
// the derived field and its sampled reference validation together.
import {readFileSync,writeFileSync,mkdirSync} from 'node:fs';
import {createHash} from 'node:crypto';
for(const p of ['engine/math.js','anatomy/shoulder-rig-data.js','engine/joint-frames.js','engine/pose-field.js','engine/shoulder-kinematics.js'])await import('../../'+p);
const input=process.argv[2]||'test-results/v076-rig-current/rig-calibration.json',fixture=JSON.parse(readFileSync(input));
const scale=fixture.scale,map=p=>[-(p[0]+14)*scale,1.7+(p[2]-1341)*scale,(p[1]+72)*scale+.08],landmarks=Object.fromEntries(Object.entries(ShoulderRigData.landmarks).map(([k,p])=>[k,map(p)]));
const model=new ShoulderKinematics.ShoulderModel({...fixture,landmarks,solverIterations:24}),elevationStep=3,planeStep=5;
const f={version:'0.76.0',method:ShoulderPoseField.method,signature:model.fieldSignature,elevationStep,planeStep,elevationMax:180,planeMax:90,elevationCount:61,planeCount:19,sourceAtlas:ShoulderRigData.atlasCommit,license:ShoulderRigData.license,referenceIterations:24,values:[]};
const time=performance.now();
for(let e=0;e<=180;e+=elevationStep){for(let p=0;p<=90;p+=planeStep)f.values.push(...model.pose({elev:e,plane:p,backend:'reference'}).parameters.map(x=>Math.fround(x)));if(e%30===0)console.log('Generated',e)}
// Conservative rigid-point error at each validation posture. Radii use the
// actual source bounding boxes, including corners outside the bone surface.
const bounds={clavicle:[[-150.03300476,-146.76600647,1331.76000977],[-11.44390011,-47.22389984,1366.85998535]],scapula:[[-168.68499756,-99.07630157,1201.5],[-59.31779861,11.35169983,1359.84997559]]};
const radius=(name,pivot)=>{const [lo,hi]=bounds[name];let r=0;for(const x of [lo[0],hi[0]])for(const y of [lo[1],hi[1]])for(const z of [lo[2],hi[2]])r=Math.max(r,Math.hypot(x-pivot[0],y-pivot[1],z-pivot[2]));return r};
const rc=radius('clavicle',ShoulderRigData.landmarks.SC),rs=radius('scapula',ShoulderRigData.landmarks.AC),M=AnatomyMath;
const angle=(a,b)=>2*Math.acos(Math.max(-1,Math.min(1,Math.abs(a.reduce((s,x,i)=>s+x*b[i],0)))));
let max=0,worst=null,count=0;
for(let e=0;e<180;e+=elevationStep){for(let p=0;p<90;p+=planeStep)for(const a of [.25,.5,.75])for(const b of [.25,.5,.75]){
 const elev=e+a*elevationStep,plane=p+b*planeStep,reference=model.pose({elev,plane,backend:'reference'}),approx=model.forward(ShoulderPoseField.sample(f,elev,plane));
 const ac=M.len(M.sub(reference.AC,approx.AC))/scale,gh=M.len(M.sub(reference.GH,approx.GH))/scale;
 const error=Math.max(2*rc*Math.sin(angle(reference.qC,approx.qC)/2),ac+2*rs*Math.sin(angle(reference.qS,approx.qS)/2),gh);
 if(error>max){max=error;worst={elev,plane,acMm:ac,ghMm:gh,rigidPointBoundMm:error}}count++;
 }if(e%30===0)console.log('Validated',e,'max mm',max)}
if(max>=1)throw Error('Pose field exceeds 1mm reference tolerance: '+JSON.stringify({max,worst}));
f.validation={kind:'sampled conservative bound for rigid points within source boxes; not clinical accuracy',samples:count,maxRigidErrorBoundMm:max,worst,totalMs:performance.now()-time,physicalPhoneTested:false};
f.referenceSourceHash=createHash('sha256').update(readFileSync('engine/shoulder-kinematics.js')).digest('hex');
const {values,...meta}=f;
// Keep the numerical payload typed at runtime without duplicating arrays.
const output='/* '+f.license+'; generated, see scripts/calibration/build-shoulder-field.mjs. */\n(function(root){"use strict";const data='+JSON.stringify(meta,null,2)+';data.values=new Float32Array('+JSON.stringify(values)+');root.ShoulderPoseFieldData=Object.freeze(data);})(globalThis);\n';
writeFileSync('anatomy/shoulder-pose-field-v076.js',output);
mkdirSync('docs/validation/v076',{recursive:true});writeFileSync('docs/validation/v076/pose-field-validation.json',JSON.stringify(meta,null,2)+'\n');
writeFileSync('tests/fixtures/shoulder-v076-field-calibration.json',JSON.stringify({...fixture,landmarks},null,2)+'\n');
console.log(f.validation);
