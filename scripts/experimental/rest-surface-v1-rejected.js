/* Rest-mesh co-rotation with a differential surface-area correction. This
 * retains the atlas surface and wide attachment masks. It never divides by a
 * fitted guide tangent/width, nor aligns a residual to a possibly folded path.
 * Passive geometric approximation, not a constitutive muscle/force model.
 */
(function(root){
 'use strict';
 const M=root.AnatomyMath,{add,sub,mul,dot,cross,len,norm,qrot,apply,mix}=M;
 function center(controls,l,t){
  l=M.clamp(l,0,controls.length-1);const i=Math.floor(l),u=l-i,v=1-u;
  const W=[v*v*v/6,(3*u*u*u-6*u*u+4)/6,(-3*u*u*u+3*u*u+3*u+1)/6,u*u*u/6],D=[-v*v/2,(9*u*u-12*u)/6,(-9*u*u+6*u+3)/6,u*u/2];
  let C=[0,0,0],Ct=[0,0,0],Cl=[0,0,0],A=[0,0,0],B=[0,0,0];
  for(let k=0;k<4;k++){const c=controls[M.clamp(i+k-1,0,controls.length-1)],p=root.FiberPath.curve(c,t);C=add(C,mul(p,W[k]));Ct=add(Ct,mul(root.FiberPath.tangent(c,t),W[k]));Cl=add(Cl,mul(p,D[k]));A=add(A,mul(c[0],W[k]));B=add(B,mul(c[3],W[k]));}
  return {C,Ct,Cl,T:norm(sub(B,A))};
 }
 function basis(B,l,t){
  const c=center(B.field.controls,l,t),T=c.T;
  let U=sub(B.reference,mul(T,dot(T,B.reference)));
  if(len(U)<1e-7){const a=Math.abs(T[0])<.8?[1,0,0]:[0,1,0];U=sub(a,mul(T,dot(T,a)));}
  U=norm(U);return {...c,U,N:norm(cross(T,U))};
 }
 function bind(B,graph){
  const n=B.weights.length,gradT=new Float32Array(n*3),gradU=new Float32Array(n*3),{offsets,indices}=graph;
  let rankDeficient=0,maxCondition=0;
  for(let i=0;i<n;i++){
   const b=basis(B,B.lane[i],B.material[i]),p=B.rest.subarray(i*3,i*3+3),rhs=[[0,0],[0,0],[0,0]];let xx=0,xy=0,yy=0;
   for(let k=offsets[i];k<offsets[i+1];k++){const j=indices[k],d=sub(B.rest.subarray(j*3,j*3+3),p),x=dot(d,b.T),y=dot(d,b.U),a=1/Math.max(dot(d,d),1e-12);xx+=a*x*x;xy+=a*x*y;yy+=a*y*y;const df=[B.weights[j]-B.weights[i],B.lane[j]-B.lane[i],B.material[j]-B.material[i]];for(let h=0;h<3;h++){rhs[h][0]+=a*x*df[h];rhs[h][1]+=a*y*df[h];}}
   const trace=xx+yy,disc=Math.hypot(xx-yy,2*xy),low=(trace-disc)/2,high=(trace+disc)/2;
   maxCondition=Math.max(maxCondition,high/Math.max(low,1e-12));if(low<1e-8*trace)rankDeficient++;
   // Tikhonov scale follows the dimensionless local covariance. No muscle-
   // specific tuning or pose sampling enters the rest differential fit.
   const ridge=Math.max(trace*1e-8,1e-12);xx+=ridge;yy+=ridge;const det=xx*yy-xy*xy;
   for(let h=0;h<3;h++){gradT[i*3+h]=(yy*rhs[h][0]-xy*rhs[h][1])/det;gradU[i*3+h]=(xx*rhs[h][1]-xy*rhs[h][0])/det;}
  }
  return {gradT,gradU,rankDeficient,maxCondition,method:'rest-surface-differential-v1'};
 }
 function rotationDerivative(q,dq,v){const qv=q.slice(0,3),dv=dq.slice(0,3);return mul(add(add(mul(cross(qv,v),dq[3]),mul(cross(dv,v),q[3])),add(cross(dv,cross(qv,v)),cross(qv,cross(dv,v)))),2);}
 function evaluate(state,p,w,l,t,index,normal=null){
  const B=state.binding,a=state.originAt?state.originAt(p):state.a,b=state.b,base=mix(apply(a,p),apply(b,p),w);
  if(w<=0||w>=1)return {point:base,normal:normal?qrot(w>=1?b.q:a.q,normal):null,stretch:1,determinant:1,bounded:false};
  if(index<0||!B.surface)throw Error('Rest-surface vertex differential missing');
  const f=basis(B,l,t),off=sub(p,f.C),raw=a.q.map((v,i)=>v*(1-w)+b.q[i]*w),ql=Math.hypot(...raw);
  if(ql<1e-6)return {point:root.ContactConstraint.project(p,base,w,state.contact),normal,invalidFrame:true,reason:'ambiguous-joint-quaternion-lift',stretch:null,determinant:null};
  const q=raw.map(x=>x/ql),dq0=b.q.map((x,i)=>x-a.q[i]),qd=q.reduce((s,x,i)=>s+x*dq0[i],0),dq=dq0.map((x,i)=>(x-q[i]*qd)/ql),A=apply(a,f.C),C=mix(A,apply(b,f.C),w),dw=add(sub(apply(b,f.C),A),rotationDerivative(q,dq,off));
  const F=(v,g)=>{const dC=add(mul(f.Cl,g[1]),mul(f.Ct,g[2]));return add(add(qrot(q,v),sub(mix(qrot(a.q,dC),qrot(b.q,dC),w),qrot(q,dC))),mul(dw,g[0]));},g=B.surface,FT=F(f.T,g.gradT.subarray(index*3,index*3+3)),FU=F(f.U,g.gradU.subarray(index*3,index*3+3)),QN=qrot(q,f.N),J=dot(cross(FT,FU),QN),gate=M.smooth(w/.12)*M.smooth((1-w)/.12),tendon=root.TendonDeformation.fraction(w,B.definition.tissues);
  let U=1,N=1,bounded=false;
  if(J>1e-8){const u0=B.sheet?1:1/Math.sqrt(J),n0=B.sheet?1/J:u0;U=M.clamp(u0,.4,1.6);N=M.clamp(n0,.4,1.6);bounded=U!==u0||N!==n0;U=1+gate*(1-tendon)*(U-1);N=1+gate*(1-tendon)*(N-1);}
  const scaled=add(off,add(mul(f.U,dot(off,f.U)*(U-1)),mul(f.N,dot(off,f.N)*(N-1))));
  const candidate=add(C,qrot(q,scaled)),point=root.ContactConstraint.project(p,candidate,w,state.contact);
  return {point,normal:normal?qrot(q,normal):null,stretch:len(FT),determinant:J*U*N,bounded,invalidFrame:J<=1e-8,areaJacobianApprox:J,volumeApproximation:'surface differential; excludes scale/contact derivatives'};
 }
 root.RestSurface=Object.freeze({center,basis,bind,evaluate});
})(globalThis);
