// Export final rigid transforms for an independent source-triangle audit.
import {readFileSync,writeFileSync,mkdirSync}from'node:fs';
for(const f of ['engine/math','anatomy/shoulder-rig-data','engine/joint-frames','anatomy/thorax-contact-data-v076','engine/convex-contact','engine/shoulder-kinematics'])await import('../../'+f+'.js');
if(process.env.QA_BONE_FIELD==='1'){await import('../../engine/pose-field.js');await import('../../anatomy/shoulder-pose-field-v076.js');}
const calibration=JSON.parse(readFileSync('docs/validation/v076/contact-rms/rig-calibration.json')),rig=new ShoulderKinematics.ShoulderModel(calibration),out=process.argv[2]||'test-results/v076-bone-surfaces';mkdirSync(out,{recursive:true});
const poses=[];for(const plane of [0,40,90])for(const elev of (process.env.QA_BONE_SWEEP==='1'?Array.from({length:181},(_,i)=>i):[0,30,60,90,120,150,180])){const p=rig.pose({elev,plane});poses.push({elev,plane,bones:p.bones,proxy:p.diagnostics})}
writeFileSync(out+'/poses.json',JSON.stringify({calibration,atlasSC:ShoulderRigData.landmarks.SC,poses},null,2));
