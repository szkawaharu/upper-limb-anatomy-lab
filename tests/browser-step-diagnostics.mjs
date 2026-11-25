// Diagnostics only: no change to the production deformation or sweep sample policy.
import assert from "node:assert/strict";
for (const name of ["math","fiber-path","volume-constraint","contact-constraint","tendon-deformation","muscle-deformation"]) await import("../engine/"+name+".js");
const M=globalThis.AnatomyMath,P=globalThis.FiberPath,E=globalThis.MuscleDeformation;
const scale=.00555;
const controls=(flat)=>Array.from({length:flat.length/12},(_,l)=>Array.from({length:4},(_,j)=>flat.slice(l*12+j*3,l*12+j*3+3)));
const dist=(a,b)=>Math.hypot(...a.map((x,i)=>x-b[i]))/scale;
export function describeStep(row,i){
 const s=row.state,x=row.samples[i],a=x.originAnchor||s.a,b=s.b,rest=s.restControls||controls(s.rest),posed=s.posedControls||controls(s.now);
 const R=P.frame(rest,x.l,x.materialT,s.reference),raw=P.frame(posed,x.l,x.materialT,M.qrot(a.q,s.reference));
 const base=M.mix(M.apply(a,x.p),M.apply(b,x.p),x.t);
 const materialBlend=s.blendRange[0]+(s.blendRange[1]-s.blendRange[0])*x.materialT;
 const binding={field:{controls:rest},reference:s.reference,blendRange:s.blendRange,sheet:s.sheet,definition:{id:row.id,tissues:{tendonStart:s.tendonStart}}};
 const state={binding,a,b,posed,contact:null};
 const pre=E.evaluate(state,x.p,x.t,x.l,null,x.materialT),post=E.evaluate({...state,contact:s.contact},x.p,x.t,x.l,null,x.materialT);
 const S=R&&raw?P.transported(R,raw,a,b,materialBlend):null;
 const refT=R?M.qrot(M.qmix(a.q,b.q,M.smooth(materialBlend)),R.t):null;
 const refU=R?M.qrot(M.qmix(a.q,b.q,M.smooth(materialBlend)),R.u):null;
 const rawLongitudinal=R&&raw?raw.tl/R.tl:null,rawTransverse=R&&raw?raw.ul/Math.max(R.ul,1e-7):null;
 const tendon=globalThis.TendonDeformation.fraction(x.t,binding.definition.tissues);
 const scales=R&&raw?globalThis.VolumeConstraint.scales(rawLongitudinal,rawTransverse,s.sheet,tendon):null;
 return {vertexId:x.vertexId,restPoint:x.p,weight:x.t,lane:x.l,materialT:x.materialT,point:x.q,originAnchor:a,insertionAnchor:b,restFrame:R,rawPosedFrame:raw,transportedFrame:S,referenceT:refT,referenceU:refU,tangentDot:refT&&raw?M.dot(refT,raw.t):null,rawLongitudinal,rawTransverse,scales,frameConfidence:S?.frameConfidence??null,effectiveTangentDot:S?.effectiveTangentDot??null,gate:M.smooth(x.t/.12)*M.smooth((1-x.t)/.12),base,preContact:pre.point,postContact:post.point,contactDeltaMm:dist(pre.point,post.point),reconstructedErrorMm:dist(post.point,x.q),invalidFrame:!!pre.invalidFrame};
}
export function trackStep(maxima,before,after,i,pose){
 const x=before.samples[i],y=after.samples[i];
 assert.equal(before.id,after.id);assert.equal(x.vertexId,y.vertexId);
 const mm=dist(x.q,y.q);
 if(!maxima[after.id]||mm>maxima[after.id].mm){
  const a=describeStep(before,i),b=describeStep(after,i);
  maxima[after.id]={id:after.id,name:after.name||after.id,side:after.side||"right",vertexId:y.vertexId,mm,...pose,before:a,after:b,stageDeltaMm:{base:dist(a.base,b.base),preContact:dist(a.preContact,b.preContact),postContact:dist(a.postContact,b.postContact)},tangentFrameDot:a.transportedFrame&&b.transportedFrame?M.dot(a.transportedFrame.u,b.transportedFrame.u):null};
 }
 return mm;
}
