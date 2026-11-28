/* Trace source physical sizes and stage-by-stage deformation of worst
 * mesh elements. Diagnostics only; no shape or contact limits changed. */
(function(root){
'use strict';
const M=root.AnatomyMath,E=root.MuscleDeformation,P=root.FiberPath;
function trace(state,report){
 const B=state.binding,scale=B.definition.source.worldUnitsPerMm;
 const ids=[...new Set([...(report.worstEdge?.vertices||[]),...(report.worstArea?.vertices||[])])];
 const vertices=ids.map(i=>{
  const p=Array.from(B.rest.slice(3*i,3*i+3)),w=B.weights[i],lane=B.lane[i],t=B.material[i],a=state.originAt?state.originAt(p):state.a;
  const base=M.mix(M.apply(a,p),M.apply(state.b,p),w),noContact=E.evaluate({...state,contact:null},p,w,lane,null,t),final=E.evaluate(state,p,w,lane,null,t);
  const R=P.frame(B.field.controls,lane,t,B.reference),S=P.frame(state.posed,lane,t,M.qrot(a.q,B.reference));
  return {i,rest:p,weight:w,lane,materialT:t,base,noContact:noContact.point,final:final.point,
   restWidthMm:R?R.ul/scale:null,restTangentMm:R?R.tl/scale:null,restOffsetMm:R?M.len(M.sub(p,R.c))/scale:null,
   posedWidthMm:S?S.ul/scale:null,posedTangentMm:S?S.tl/scale:null,stretch:final.stretch,bounded:final.bounded,
   originBoneTransform:a};
 });
 const lookup=new Map(vertices.map(v=>[v.i,v]));
 const edge=(ids,key)=>M.len(M.sub(lookup.get(ids[0])[key],lookup.get(ids[1])[key]))/scale;
 const area=(ids,key)=>{const p=ids.map(i=>lookup.get(i)[key]);return M.len(M.cross(M.sub(p[1],p[0]),M.sub(p[2],p[0])))/(2*scale*scale)};
 const stages=['rest','base','noContact','final'];
 return {id:report.id,fieldMethod:B.fieldMethod,vertices,worstEdge:report.worstEdge?{...report.worstEdge,stageLengthsMm:Object.fromEntries(stages.map(k=>[k,edge(report.worstEdge.vertices,k)]))}:null,worstArea:report.worstArea?{...report.worstArea,stageAreasMm2:Object.fromEntries(stages.map(k=>[k,area(report.worstArea.vertices,k)]))}:null};
}
root.MuscleOutlierTrace=Object.freeze({trace});
})(globalThis);
