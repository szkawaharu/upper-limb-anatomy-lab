/* Experimental shared rest-surface transport and a scalar material-volume
 * projection. All shape information comes from the source surface and its
 * static rest data. The scalar is solved from a closed-source surface integral,
 * never fitted by muscle name/pose or interpreted as activation/force.
 */
(function(root){
 'use strict';
 const M=root.AnatomyMath,{add,sub,mul,mix,qrot,apply,dot,cross}=M;
 function raw(state,index){
  const B=state.binding,p=B.rest.subarray(index*3,index*3+3),w=B.weights[index],a=state.originAt?state.originAt(p):state.a,b=state.b,C=B.transport.centers.subarray(index*3,index*3+3);
  if(w===0||w===1)return {p:apply(w?b:a,p),d:[0,0,0],q:(w?b:a).q};
  const q=M.qblendLift(a.q,b.q,w),off=sub(p,C),center=mix(apply(a,C),apply(b,C),w),gate=M.smooth(w/.12)*M.smooth((1-w)/.12)*(1-root.TendonDeformation.fraction(w,B.definition.tissues));
  return {p:add(center,qrot(q,off)),d:mul(qrot(q,B.transport.offsets.subarray(index*3,index*3+3)),gate),q};
 }
 function project(state,index,mu){const B=state.binding,r=raw(state,index),p=add(r.p,mul(r.d,mu));return root.ContactConstraint.project(B.rest.subarray(index*3,index*3+3),p,B.weights[index],state.contact);}
 function volume(state){
  const B=state.binding,Q=B.transport.quadrature;
  if(!Q)return {mu:0,status:'unavailable',reason:'source-volume-not-evaluable',iterations:0,ratio:null};
  if(state.liftSingular)return {mu:0,status:'unavailable',reason:'ambiguous-joint-quaternion-lift',iterations:0,ratio:null};
  const ids=Q.sourceVertices,o=apply(state.a,Q.sourceCenter),raws=ids.map(i=>raw(state,i)),points=new Float64Array(ids.length*3),reference=Q.referenceVolume;
  const insertionDelta=dot(sub(apply(state.b,Q.sourceCenter),o),qrot(state.b.q,Q.pureBVector));
  const value=mu=>{
   for(let j=0;j<ids.length;j++){const i=ids[j],r=raws[j],p=add(r.p,mul(r.d,mu)),q=root.ContactConstraint.project(B.rest.subarray(i*3,i*3+3),p,B.weights[i],state.contact);for(let k=0;k<3;k++)points[j*3+k]=q[k]-o[k];}
   let delta=insertionDelta;for(let k=0;k<Q.multipliers.length;k++){const i=k*3,a=points.subarray(Q.faces[i]*3,Q.faces[i]*3+3),b=points.subarray(Q.faces[i+1]*3,Q.faces[i+1]*3+3),c=points.subarray(Q.faces[i+2]*3,Q.faces[i+2]*3+3);delta+=Q.multipliers[k]*(dot(a,cross(b,c))-Q.controlDet[k]);}return reference+delta/6;
  };
  let lo=-.6,hi=.6,flo=value(lo)-reference,fhi=value(hi)-reference,f0=value(0)-reference;
  if(flo>0||fhi<0||f0<flo||f0>fhi)return {mu:0,status:'unavailable',reason:'volume-projection-unbracketed',iterations:0,ratio:1+f0/reference,bounds:[1+flo/reference,1+fhi/reference]};
  if(Math.abs(f0)<reference*1e-6)return {mu:0,status:'solved',iterations:0,ratio:1+f0/reference};
  let mu=0,f=f0;
  for(let i=0;i<16;i++){
   mu=(lo*fhi-hi*flo)/(fhi-flo);f=value(mu)-reference;
   if(Math.abs(f)<reference*1e-6)return {mu,status:'solved',iterations:i+1,ratio:1+f/reference};
   if(f<0){lo=mu;flo=f;}else{hi=mu;fhi=f;}
  }
  return {mu,status:'iteration-limit',iterations:16,ratio:1+f/reference};
 }
 function deform(state,mu=state.materialVolume?.mu||0){const B=state.binding,out=new Float32Array(B.rest.length);for(let i=0;i<B.weights.length;i++)out.set(project(state,i,mu),i*3);return out;}
 root.RestTransport=Object.freeze({raw,project,volume,deform});
})(globalThis);
