// Build an OFFLINE candidate field for the screened continuous blend.
// Existing product field remains unchanged. Reference is the same
// conservative-field + independent loft reference blend as the source audit.
import {readFileSync,writeFileSync,mkdirSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {readBlend} from './continuous-shoulder-blend.mjs';
for(const p of ['engine/math','anatomy/shoulder-rig-data','engine/joint-frames','engine/pose-field','anatomy/shoulder-pose-field-v076','anatomy/thorax-contact-data-v076','engine/convex-contact','engine/shoulder-kinematics'])await import('../../'+p+'.js');
const fixture=JSON.parse(readFileSync('docs/validation/v076/contact-rms/rig-calibration.json')),M=AnatomyMath;
const current=new ShoulderKinematics.ShoulderModel(fixture),local=new ShoulderKinematics.ShoulderModel(fixture);
if(!current.field)throw Error('Current reference field signature mismatch');
local.body=null;local.field=null;
const blend=readBlend(),base=current.field;
// Low-elevation guard derived from the observed worst worsening at 35 degrees
// and the last worsening at 40 degrees. It preserves the established neutral
// behavior. It is an experiment, not a physiological transition angle.
const gate=e=>{const t=Math.max(0,Math.min(1,(e-40)/5));return t*t*(3-2*t)};
const coefficient=(e,p)=>blend.sample(e,p);
const refCache=new Map();
function reference(e,p){
 const k=e+':'+p;if(refCache.has(k))return refCache.get(k);
 const a=current.pose({elev:e,plane:p}).parameters;
 const t=coefficient(e,p),b=t?local.pose({elev:e,plane:p,backend:'reference'}).parameters:a;
 const v=a.map((x,j)=>t*(b[j]-x));refCache.set(k,v);return v;
}
let E=[...new Set([...base.elevationKnots,...blend.nodes,40,45])].sort((a,b)=>a-b);
let P=[...new Set([...base.planeKnots,...blend.planes])].sort((a,b)=>a-b);
const bounds={clavicle:[[-150.03300476,-146.76600647,1331.76000977],[-11.44390011,-47.22389984,1366.85998535]],scapula:[[-168.68499756,-99.07630157,1201.5],[-59.31779861,11.35169983,1359.84997559]]};
const radius=(name,pivot)=>{const [l,h]=bounds[name];return Math.max(...[l[0],h[0]].flatMap(x=>[l[1],h[1]].flatMap(y=>[l[2],h[2]].map(z=>Math.hypot(x-pivot[0],y-pivot[1],z-pivot[2])))))};
const rc=radius('clavicle',ShoulderRigData.landmarks.SC),rs=radius('scapula',ShoulderRigData.landmarks.AC);
const qdist=(a,b)=>2*Math.sqrt(Math.max(0,1-Math.min(1,Math.abs(a.reduce((s,x,i)=>s+x*b[i],0)))**2));
function error(a,b){
 const x=current.forward(a),y=current.forward(b),ac=M.len(M.sub(x.AC,y.AC))/fixture.scale,gh=M.len(M.sub(x.GH,y.GH))/fixture.scale;
 return Math.max(rc*qdist(x.qC,y.qC),ac+rs*qdist(x.qS,y.qS),gh);
}
const levels=[];let f,max,worst;
for(let level=0;level<4;level++){
 f={elevationKnots:E,planeKnots:P,values:E.flatMap(e=>P.flatMap(p=>reference(e,p).map(Math.fround)))};
 max=0;worst=null;let samples=0;const badE=new Set(),badP=new Set();
 for(let i=0;i<E.length-1;i++)for(let j=0;j<P.length-1;j++)for(const u of [.25,.5,.75])for(const v of [.25,.5,.75]){
  const e=E[i]+u*(E[i+1]-E[i]),p=P[j]+v*(P[j+1]-P[j]);
  const a=current.pose({elev:e,plane:p}).parameters,g=gate(e),d=reference(e,p),z=ShoulderPoseField.sample(f,e,p);
  const err=error(a.map((x,k)=>x+g*d[k]),a.map((x,k)=>x+g*z[k]));
  if(err>max){max=err;worst={elev:e,plane:p,rigidPointBoundMm:err}}
  if(err>=1){badE.add((E[i]+E[i+1])/2);badP.add((P[j]+P[j+1])/2)}
  samples++;
 }
 levels.push({level,nodes:E.length*P.length,samples,max,worst});
 console.log('PREVIEW_FIELD_LEVEL '+JSON.stringify(levels.at(-1)));
 if(max<1)break;
 E=[...new Set([...E,...badE])].sort((a,b)=>a-b);P=[...new Set([...P,...badP])].sort((a,b)=>a-b);
}
if(!(max<1))throw Error('Candidate interpolation fails unchanged 1mm engineering tolerance');
const deps=['engine/math.js','anatomy/shoulder-rig-data.js','engine/joint-frames.js','engine/pose-field.js','anatomy/shoulder-pose-field-v076.js','anatomy/thorax-contact-data-v076.js','engine/convex-contact.js','engine/shoulder-kinematics.js','scripts/calibration/continuous-shoulder-blend.mjs','scripts/calibration/build-shoulder-blend-preview.mjs','docs/validation/v076/shoulder-range-screening.json','docs/validation/v076/contact-rms/rig-calibration.json'];
const out=process.env.QA_PREVIEW_OUT||'test-results/v076-blend-preview';mkdirSync(out,{recursive:true});
const meta={kind:'offline additive six-parameter correction preview; not production',representation:'current conservative parameters + smoothstep(40,45,elev) * interpolated delta',baseSignature:current.fieldSignature,lowElevationGate:[40,45],sourceRun:blend.sourceRun,sourceHashes:Object.fromEntries(deps.map(p=>[p,createHash('sha256').update(readFileSync(p)).digest('hex')])),levels,maxRigidErrorBoundMm:max,worst,productAccepted:false,anatomicalAcceptance:false};
writeFileSync(out+'/field.json',JSON.stringify({...meta,...f})+'\n');
writeFileSync(out+'/summary.json',JSON.stringify(meta,null,2)+'\n');
console.log('PREVIEW_FIELD_SUMMARY '+JSON.stringify(meta));
