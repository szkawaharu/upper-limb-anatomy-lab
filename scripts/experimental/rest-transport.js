/* Experimental shared rest-surface transport and a scalar material-volume
 * projection. All shape information comes from the source surface and its
 * static rest data. The scalar is solved from a closed-source surface integral,
 * never fitted by muscle name/pose or interpreted as activation/force.
 */
(function(root){
 'use strict';
 const M=root.AnatomyMath,{add,sub,mul,mix,qrot,apply,dot,cross}=M;
 function surfaceRaw(state,index){
  const B=state.binding,p=B.rest.subarray(index*3,index*3+3),w=B.weights[index],a=state.originAt?state.originAt(p):state.a,b=state.b,C=B.transport.centers.subarray(index*3,index*3+3);
  if(w===0||w===1)return {p:apply(w?b:a,p),d:[0,0,0],q:(w?b:a).q};
  const q=M.qblendLift(a.q,b.q,w),off=sub(p,C),center=mix(apply(a,C),apply(b,C),w),gate=M.smooth(w/.12)*M.smooth((1-w)/.12)*(1-root.TendonDeformation.fraction(w,B.definition.tissues));
  return {p:add(center,qrot(q,off)),d:mul(qrot(q,B.transport.offsets.subarray(index*3,index*3+3)),gate),q};
 }
 function chordRaw(state,index){
  const B=state.binding,p=B.rest.subarray(index*3,index*3+3),w=B.weights[index],a=state.originAt?state.originAt(p):state.a,b=state.b;
  if(w===0||w===1)return {p:apply(w?b:a,p),d:[0,0,0],q:(w?b:a).q};
  const u=M.clamp((w-B.blendRange[0])/(B.blendRange[1]-B.blendRange[0]),0,1),c=root.FiberPath.laneCurve(B.field.controls,B.lane[index]),O=c[0],I=c[3],v=sub(I,O),T=M.norm(v),op=apply(a,O),ip=apply(b,I),V=sub(ip,op),S=M.norm(V),q=M.qblendLift(a.q,b.q,w),ref=B.reference;
  let U=sub(ref,mul(T,dot(T,ref)));if(M.len(U)<1e-7){const e=Math.abs(T[0])<.8?[1,0,0]:[0,1,0];U=sub(e,mul(T,dot(T,e)));}U=M.norm(U);const N=M.norm(cross(T,U)),controls=B.field.controls,L=controls.at(-1),F=controls[0],um=mix(sub(L[0],F[0]),sub(L[3],F[3]),u),up=mix(qrot(a.q,sub(L[0],F[0])),qrot(b.q,sub(L[3],F[3])),u),rwidth=M.len(sub(um,mul(T,dot(um,T)))),pwidth=M.len(sub(up,mul(S,dot(up,S)))),sc=root.VolumeConstraint.scales(M.len(V)/Math.max(M.len(v),1e-9),rwidth>1e-7?pwidth/rwidth:1,B.sheet,root.TendonDeformation.fraction(w,B.definition.tissues)),off=sub(p,mix(O,I,u)),scaled=add(mul(T,dot(off,T)*sc.t),add(mul(U,dot(off,U)*sc.u),mul(N,dot(off,N)*sc.n))),QT=qrot(q,T),QU=qrot(q,U),rot=x=>M.align(QT,S,qrot(q,x),QU),candidate=add(mix(op,ip,u),rot(scaled)),gate=M.smooth(w/.12)*M.smooth((1-w)/.12),base=mix(apply(a,p),apply(b,p),w),d=mul(rot(B.transport.offsets.subarray(index*3,index*3+3)),gate*(1-root.TendonDeformation.fraction(w,B.definition.tissues)));
  return {p:mix(base,candidate,gate),d,q};
 }
 function dualRaw(state,index){
  const B=state.binding,p=B.rest.subarray(index*3,index*3+3),w=B.weights[index],a=state.originAt?state.originAt(p):state.a,b=state.b;
  if(w===0||w===1)return {p:apply(w?b:a,p),d:[0,0,0],q:(w?b:a).q};
  const dual=A=>M.qmul([...sub(A.p,qrot(A.q,A.rest)),0],A.q).map(x=>x*.5),da=dual(a),db=dual(b),real=a.q.map((x,i)=>x*(1-w)+b.q[i]*w),n=Math.hypot(...real);
  if(n<1e-8)throw Error('Ambiguous continuous dual-quaternion lift');
  const q=real.map(x=>x/n),d=da.map((x,i)=>(x*(1-w)+db[i]*w)/n),v=M.qmul(d,[-q[0],-q[1],-q[2],q[3]]),gate=M.smooth(w/.12)*M.smooth((1-w)/.12)*(1-root.TendonDeformation.fraction(w,B.definition.tissues));
  return {p:add(qrot(q,p),v.slice(0,3).map(x=>2*x)),d:mul(qrot(q,B.transport.offsets.subarray(index*3,index*3+3)),gate),q};
 }
 const raw=(state,index)=>state.binding.transport.kind==='chord'?chordRaw(state,index):state.binding.transport.kind==='dqs'?dualRaw(state,index):surfaceRaw(state,index);
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
