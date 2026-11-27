// Offline candidate sweep for the visible controls in the 61.33 s recording.
 // Reduce excessive separation caused by the whole-thorax convex enclosure
 // only when exact source-triangle checks subsequently permit it.
import {readFileSync,writeFileSync,mkdirSync} from 'node:fs';
for(const f of ['engine/math','anatomy/shoulder-rig-data','engine/joint-frames','engine/pose-field','anatomy/shoulder-pose-field-v076','anatomy/thorax-contact-data-v076','engine/convex-contact','engine/shoulder-kinematics','engine/rig-diagnostics'])await import('../../'+f+'.js');
const calibration=JSON.parse(readFileSync('docs/validation/v076/contact-rms/rig-calibration.json'));
const current=new ShoulderKinematics.ShoulderModel(calibration),local=new ShoulderKinematics.ShoulderModel(calibration);
local.body=null;local.field=null; // comparison endpoint ONLY
const rows=[];
const sweep=process.env.QA_SHOULDER_PLANE!==undefined;
const plane=Number(process.env.QA_SHOULDER_PLANE),step=Number(process.env.QA_SHOULDER_STEP||15);
if(sweep&&(!Number.isFinite(plane)||plane<0||plane>90||!Number.isInteger(step)||step<1||180%step))throw Error('Invalid sweep');
const inputs=sweep?Array.from({length:180/step+1},(_,i)=>({elev:i*step,plane})):
 [{elev:0,plane:0},{elev:70,plane:39},{elev:180,plane:44},{elev:150,plane:90}];
for(const input of inputs){
 const a=current.pose(input),b=local.pose({...input,backend:'reference'});
 for(const alpha of input.elev===0?[0]:sweep?[0,.125,.25]:[0,.25,.5,.75,1]){
  const params=a.parameters.map((x,j)=>x+alpha*(b.parameters[j]-x));
  const p=current.forward(params),L=current.landmarks;
  // Keep the linked pivots and absolute humeral orientation; no independent
  // translation or angle clipping of the scapula.
  const clavicle={q:p.qC,p:L.SC,rest:L.SC},
   scapula=JointFrames.attach(clavicle,L.AC,AnatomyMath.qmul(JointFrames.conjugate(p.qC),p.qS)),
   humerus=JointFrames.attach(scapula,L.GH,AnatomyMath.qmul(JointFrames.conjugate(p.qS),a.qH));
  const bones={clavicle,scapula,humerus};
  const connections=RigDiagnostics.connections(current,bones);
  if(connections.some(x=>x.sharedPivotErrorMm>1e-7))throw Error('Unlinked candidate');
  rows.push({...input,alpha,parameters:params,bones,connections,controlBackend:a.diagnostics.backend,
   priorDeviationL2Deg:Math.hypot(...params.map((x,j)=>x-a.prior[j])),productAccepted:false});
 }
}
const out=process.env.QA_SHOULDER_OUT||'test-results/v076-recording-shoulder';mkdirSync(out,{recursive:true});
writeFileSync(out+'/candidates.json',JSON.stringify({calibration,atlasSC:ShoulderRigData.landmarks.SC,rows,limitations:'One canonical shoulder at each observed UI setting. No individual anatomical measurements; no continuous field or muscle acceptance.'},null,2)+'\n');
console.log(JSON.stringify({candidates:rows.length,inputs,sweep,productAccepted:false}));
