/* Independent per-bone boundary targets for the experimental levator mesh.
 * A near-surface candidate is not a measured anatomical footprint.
 * This module prescribes only attachment vertices; no interior deformation,
 * no new pins, no source geometry modification. Canonical-left transforms
 * must be supplied by caller in the same chart as the rest mesh. */
(function(root){
'use strict';
const keys=['C1','C2','C3','C4','scapula'];
function registry(data){
 const n=data.rest.length/3;
 if(!Number.isInteger(n)||data.nearestC.length!==n)throw Error('Invalid source registry');
 const label=new Int8Array(n).fill(-1),counts=[0,0,0,0,0];
 for(const [ids,isInsertion]of [[data.origin,false],[data.insertion,true]])for(const i of ids){
  if(!Number.isInteger(i)||i<0||i>=n||label[i]>=0)throw Error('Conflicting attachment index');
  const k=isInsertion?4:data.nearestC[i];
  if(!Number.isInteger(k)||k<0||k>(isInsertion?4:3))throw Error('Unknown origin bone');
  label[i]=k;counts[k]++;
 }
 if(counts.some(x=>x===0))throw Error('Incomplete five-bone candidate');
 return {label,counts,keys:keys.slice(),productAccepted:false};
}
function targets(rest,reg,bones){
 const M=root.AnatomyMath,rows=[];
 for(const key of keys){const a=bones[key];if(!a||![...a.q,...a.p,...a.rest].every(Number.isFinite)||Math.abs(Math.hypot(...a.q)-1)>1e-6)throw Error('Missing or invalid rigid frame '+key)}
 for(let i=0;i<reg.label.length;i++){
  const k=reg.label[i];if(k<0)continue;
  const p=Array.from(rest.slice(3*i,3*i+3));
  rows.push({i,bone:keys[k],point:M.apply(bones[keys[k]],p)});
 }
 return rows;
}
function audit(state,data,reg,bones,now){
 const B=state.binding,scale=B.definition.source.worldUnitsPerMm;
 if(data.rest.length!==B.rest.length||data.faces.length!==B.faces.length||data.rest.some((x,i)=>x!==B.rest[i])||data.faces.some((x,i)=>x!==B.faces[i]))throw Error('Registry source mismatch');
 const rows=targets(B.rest,reg,bones),byBone=Object.fromEntries(keys.map(k=>[k,{count:0,maxLegacyOriginErrorMm:0,maxDeformedErrorMm:0,sumSquaredMm:0}]));
 for(const {i,bone,point}of rows){
  const p=Array.from(B.rest.slice(3*i,3*i+3)),a=bone==='scapula'?state.b:state.originAt?state.originAt(p):state.a;
  const targetError=root.AnatomyMath.len(root.AnatomyMath.sub(root.AnatomyMath.apply(a,p),point))/scale;
  const error=root.AnatomyMath.len(root.AnatomyMath.sub(Array.from(now.slice(3*i,3*i+3)),point))/scale;
  const q=byBone[bone];q.count++;q.maxLegacyOriginErrorMm=Math.max(q.maxLegacyOriginErrorMm,targetError);q.maxDeformedErrorMm=Math.max(q.maxDeformedErrorMm,error);q.sumSquaredMm+=error*error;
 }
 for(const q of Object.values(byBone)){q.rmsDeformedErrorMm=Math.sqrt(q.sumSquaredMm/q.count);delete q.sumSquaredMm}
 return {byBone,productAccepted:false,anatomicalAcceptance:false};
}
root.LevatorBoundaryRegistry=Object.freeze({registry,targets,audit});
})(globalThis);
