import {readFileSync,writeFileSync} from 'node:fs';
import {rebuildChart,restSurfacePoint} from './chart-candidates.mjs';
for(const f of ['engine/math','anatomy/shoulder-rig-data','engine/joint-frames','anatomy/thorax-contact-data-v076','engine/convex-contact','engine/shoulder-kinematics','engine/fiber-path','engine/volume-constraint','engine/contact-constraint','engine/tendon-deformation','engine/muscle-deformation','engine/deformation-audit','engine/material-coordinates','scripts/experimental/rest-surface'])await import('../../'+f+'.js');
const M=AnatomyMath,E=MuscleDeformation,A=DeformationAudit,bindings=JSON.parse(readFileSync('test-results/v076-source-bindings/bindings.json')),
rig=new ShoulderKinematics.ShoulderModel(JSON.parse(readFileSync('docs/validation/v076/contact-rms/rig-calibration.json')));
const output=[];
for(const kind of ['differential']){
 let maxEdge=0,maxArea=0,maxMove=0,worstEdge=null,worstMove=null;const rows=[]; const point=(s,p,w,l,t,i)=>RestSurface.evaluate(s,p,w,l,t,i).point,deform=s=>Float32Array.from(Array.from({length:s.binding.weights.length},(_,i)=>point(s,s.binding.rest.subarray(i*3,i*3+3),s.binding.weights[i],s.binding.lane[i],s.binding.material[i],i)).flat());
 for(const raw of bindings){
  const B=rebuildChart({...raw,weights:Float32Array.from(raw.weights),rest:Float32Array.from(raw.rest),faces:Uint32Array.from(raw.faces)},kind);
  const h=readFileSync('test-results/v076-source-thickness/'+B.spec.id+'.f32');B.thickness=new Float32Array(h.buffer.slice(h.byteOffset,h.byteOffset+h.byteLength));
  B.surface=RestSurface.bind(B,MaterialCoordinates.adjacency(Array.from({length:B.weights.length},(_,i)=>B.rest.subarray(i*3,i*3+3)),B.faces));
  for(const plane of [0,40,90]){
   let previous=null;
   for(let elev=0;elev<=180;elev+=3){
    const p=rig.pose({elev,plane}),bone=id=>id==='trunk'?M.identity():p.bones[id],state=E.pose(B,bone(B.spec.a),bone(B.spec.b)),ids=Array.from({length:Math.ceil(B.weights.length/Math.ceil(B.weights.length/80))},(_,i)=>i*Math.ceil(B.weights.length/80)),now=ids.flatMap(i=>point(state,B.rest.subarray(i*3,i*3+3),B.weights[i],B.lane[i],B.material[i],i));
    if(previous)for(let i=0;i<now.length;i+=3){const d=Math.hypot(now[i]-previous[i],now[i+1]-previous[i+1],now[i+2]-previous[i+2])/.00555;if(d>maxMove){maxMove=d;worstMove={id:B.spec.id,vertex:ids[i/3],plane,elev}}}
    if([0,30,60,90,120,150,180].includes(elev)){
     const r=A.audit(state,deform(state));if(r.edgeRatio[1]>maxEdge){maxEdge=r.edgeRatio[1];worstEdge={id:B.spec.id,plane,elev}}maxArea=Math.max(maxArea,r.areaRatio[1]);rows.push({id:B.spec.id,plane,elev,edge:r.edgeRatio,area:r.areaRatio,volumeRatio:r.volumeRatio,topology:r.topology,volumeStatus:r.volumeStatus});
    }
    previous=now;
   }
  }
 }
 output.push({kind,maxEdge,maxArea,maxMove,worstEdge,worstMove,rows});console.log(JSON.stringify(output.at(-1)));
}
writeFileSync('test-results/v076-material-baseline/source-normal-candidate.json',JSON.stringify(output,null,2));
