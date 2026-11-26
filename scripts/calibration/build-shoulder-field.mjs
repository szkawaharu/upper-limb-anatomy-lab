// Offline generation. Run after an actual-atlas browser export, then commit
// the derived field and its sampled reference validation together.
import {readFileSync,writeFileSync,mkdirSync} from 'node:fs';
import {createHash} from 'node:crypto';
for(const p of ['engine/math.js','anatomy/shoulder-rig-data.js','engine/joint-frames.js','engine/pose-field.js','anatomy/thorax-contact-data-v076.js','engine/convex-contact.js','engine/shoulder-kinematics.js'])await import('../../'+p);
const input=process.argv[2]||'test-results/v076-rig-current/rig-calibration.json',fixture=JSON.parse(readFileSync(input));
const scale=fixture.scale,map=p=>[-(p[0]+14)*scale,1.7+(p[2]-1341)*scale,(p[1]+72)*scale+.08],landmarks=Object.fromEntries(Object.entries(ShoulderRigData.landmarks).map(([k,p])=>[k,map(p)]));
const model=new ShoulderKinematics.ShoulderModel({...fixture,landmarks,solverIterations:48}),M=AnatomyMath;
let elevationKnots=Array.from({length:61},(_,i)=>i*3),planeKnots=Array.from({length:19},(_,i)=>i*5);
// Refine the cell from the recorded failed uniform-grid audit, then refine
// any further measured failures. Threshold remains 1mm, on rigid-point bounds.
elevationKnots.push(127.5);planeKnots.push(32.5);elevationKnots.sort((a,b)=>a-b);planeKnots.sort((a,b)=>a-b);
const f={version:'0.76.0',method:ShoulderPoseField.method,signature:model.fieldSignature,elevationMax:180,planeMax:90,sourceAtlas:ShoulderRigData.atlasCommit,license:ShoulderRigData.license,referenceIterations:48,values:[]};
const sourceHash=createHash('sha256').update(readFileSync('engine/shoulder-kinematics.js')).digest('hex'),cacheFile='test-results/v076-source-envelope-reference-cache.json';
const dependencies=['engine/math.js','anatomy/shoulder-rig-data.js','engine/joint-frames.js','engine/convex-contact.js','anatomy/thorax-contact-data-v076.js','engine/shoulder-kinematics.js'];
const dependencyHashes=Object.fromEntries(dependencies.map(p=>[p,createHash('sha256').update(readFileSync(p)).digest('hex')]));
const cacheKey=createHash('sha256').update(JSON.stringify([dependencyHashes,model.fieldSignature])).digest('hex');
let cache={};try{const saved=JSON.parse(readFileSync(cacheFile));if(saved.key===cacheKey)cache=saved.values}catch(_){}
const reference=(elev,plane)=>{const key=elev+':'+plane;if(!cache[key])cache[key]=model.pose({elev,plane,backend:'reference'}).parameters;return model.forward(cache[key]);};
const save=()=>writeFileSync(cacheFile,JSON.stringify({key:cacheKey,dependencyHashes,values:cache}));
const time=performance.now();
// Conservative rigid-point error at each validation posture. Radii use the
// actual source bounding boxes, including corners outside the bone surface.
const bounds={clavicle:[[-150.03300476,-146.76600647,1331.76000977],[-11.44390011,-47.22389984,1366.85998535]],scapula:[[-168.68499756,-99.07630157,1201.5],[-59.31779861,11.35169983,1359.84997559]]};
const radius=(name,pivot)=>{const [lo,hi]=bounds[name];let r=0;for(const x of [lo[0],hi[0]])for(const y of [lo[1],hi[1]])for(const z of [lo[2],hi[2]])r=Math.max(r,Math.hypot(x-pivot[0],y-pivot[1],z-pivot[2]));return r};
const rc=radius('clavicle',ShoulderRigData.landmarks.SC),rs=radius('scapula',ShoulderRigData.landmarks.AC);
const angle=(a,b)=>2*Math.acos(Math.max(-1,Math.min(1,Math.abs(a.reduce((s,x,i)=>s+x*b[i],0)))));
let max=0,worst=null,count=0,levels=[];
for(let level=0;level<5;level++){
 f.elevationKnots=elevationKnots;f.planeKnots=planeKnots;f.elevationCount=elevationKnots.length;f.planeCount=planeKnots.length;f.values=[];
 for(const e of elevationKnots){for(const p of planeKnots){reference(e,p);f.values.push(...cache[e+':'+p].map(Math.fround));}if(e%30===0){save();console.log('Generated',level,e);}}
 max=0;worst=null;count=0;const badE=new Set(),badP=new Set();
 for(let i=0;i<elevationKnots.length-1;i++){
  for(let j=0;j<planeKnots.length-1;j++)for(const a of [.25,.5,.75])for(const b of [.25,.5,.75]){
   const elev=elevationKnots[i]+a*(elevationKnots[i+1]-elevationKnots[i]),plane=planeKnots[j]+b*(planeKnots[j+1]-planeKnots[j]),ref=reference(elev,plane),approx=model.forward(ShoulderPoseField.sample(f,elev,plane));
   const ac=M.len(M.sub(ref.AC,approx.AC))/scale,gh=M.len(M.sub(ref.GH,approx.GH))/scale,error=Math.max(2*rc*Math.sin(angle(ref.qC,approx.qC)/2),ac+2*rs*Math.sin(angle(ref.qS,approx.qS)/2),gh);
   if(error>max){max=error;worst={elev,plane,acMm:ac,ghMm:gh,rigidPointBoundMm:error}}if(error>=1){badE.add((elevationKnots[i]+elevationKnots[i+1])/2);badP.add((planeKnots[j]+planeKnots[j+1])/2);}count++;
  }
  if(elevationKnots[i]%30===0){save();console.log('Validated',level,elevationKnots[i],'max mm',max);}
 }
 levels.push({level,nodes:f.elevationCount*f.planeCount,samples:count,max,worst});save();if(max<1)break;
 elevationKnots=[...new Set([...elevationKnots,...badE])].sort((a,b)=>a-b);planeKnots=[...new Set([...planeKnots,...badP])].sort((a,b)=>a-b);
}
if(max>=1)throw Error('Pose field exceeds unchanged 1mm reference tolerance: '+JSON.stringify({max,worst,levels}));
f.validation={kind:'sampled conservative bound for rigid points within source boxes; not clinical accuracy',samples:count,adaptiveLevels:levels,referenceCacheEntries:Object.keys(cache).length,maxRigidErrorBoundMm:max,worst,totalMs:performance.now()-time,physicalPhoneTested:false};
f.referenceSourceHash=sourceHash;f.referenceDependencyHashes=dependencyHashes;
if(createHash('sha256').update(readFileSync('engine/shoulder-kinematics.js')).digest('hex')!==sourceHash)throw Error('Reference solver changed during generation');
const {values,...meta}=f;
// Keep the numerical payload typed at runtime without duplicating arrays.
const output='/* '+f.license+'; generated, see scripts/calibration/build-shoulder-field.mjs. */\n(function(root){"use strict";const data='+JSON.stringify(meta,null,2)+';data.values=new Float32Array('+JSON.stringify(values)+');root.ShoulderPoseFieldData=Object.freeze(data);})(globalThis);\n';
writeFileSync('anatomy/shoulder-pose-field-v076.js',output);
mkdirSync('docs/validation/v076',{recursive:true});writeFileSync('docs/validation/v076/pose-field-validation.json',JSON.stringify(meta,null,2)+'\n');
writeFileSync('tests/fixtures/shoulder-v076-field-calibration.json',JSON.stringify({...fixture,landmarks},null,2)+'\n');
console.log(f.validation);
