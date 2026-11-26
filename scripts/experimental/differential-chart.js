/* Experimental ridge-regularized two-coordinate differential. The ridge is
 * derived from incident source edge spacing, never a per-muscle fitted scale.
 * This is a geometric approximation, not a volumetric elasticity solver.
 */
(function(root){
 'use strict';
 const M=root.AnatomyMath,P=root.FiberPath,E=root.MuscleDeformation,{add,sub,mul,dot,len,norm,qrot,mix,apply}=M;
 function prepare(B){
  const sum=new Float64Array(B.weights.length),count=new Uint32Array(B.weights.length);
  for(let k=0;k<B.faces.length;k+=3)for(let j=0;j<3;j++){
   const i=B.faces[k+j],n=B.faces[k+(j+1)%3],a=B.rest.subarray(i*3,i*3+3),b=B.rest.subarray(n*3,n*3+3),r=dot(sub(a,b),sub(a,b));sum[i]+=r;sum[n]+=r;count[i]++;count[n]++;
  }
  B.chartSpacing=Float32Array.from(sum,(x,i)=>count[i]?Math.sqrt(x/count[i]):0);return B;
 }
 function jacobian(controls,l,t){
  const c=P.laneCurve(controls,l),lo=M.clamp(l-.25,0,controls.length-1),hi=M.clamp(l+.25,0,controls.length-1);
  return {c:P.curve(c,t),t:P.tangent(c,t),u:mul(sub(P.curve(P.laneCurve(controls,hi),t),P.curve(P.laneCurve(controls,lo),t)),1/Math.max(hi-lo,1e-6))};
 }
 function evaluate(state,index){
  const B=state.binding,p=B.rest.subarray(index*3,index*3+3),w=B.weights[index],a=state.originAt?state.originAt(p):state.a,b=state.b,base=mix(apply(a,p),apply(b,p),w),l=B.lane[index],t=B.material[index];
  if(w===0||w===1)return {point:base,stretch:1,determinant:1,bounded:false};
  const R=P.frame(B.field.controls,l,t,B.reference),r=jacobian(B.field.controls,l,t),s=jacobian(state.posed,l,t),sigma=B.chartSpacing[index],q=M.qblendLift(a.q,b.q,M.smooth(B.blendRange[0]+(B.blendRange[1]-B.blendRange[0])*t));
  if(!R||state.liftSingular)return {point:base,invalidFrame:true,stretch:null,determinant:null,reason:'invalid-rest-chart-or-quaternion-lift'};
  const g00=dot(r.t,r.t)+sigma*sigma,g11=dot(r.u,r.u)+sigma*sigma,g01=dot(r.t,r.u),det=g00*g11-g01*g01,K0=sub(s.t,qrot(q,r.t)),K1=sub(s.u,qrot(q,r.u));
  if(!(det>0))return {point:base,invalidFrame:true,stretch:null,determinant:null,reason:'singular-regularized-chart'};
  const differential=v=>{const d0=dot(r.t,v),d1=dot(r.u,v),z0=(g11*d0-g01*d1)/det,z1=(g00*d1-g01*d0)/det;return add(qrot(q,v),add(mul(K0,z0),mul(K1,z1)));},FT=differential(R.t),FU=differential(R.u),L=len(FT),T=norm(FT),shear=dot(FU,T),perp=sub(FU,mul(T,shear)),width=len(perp),U=norm(perp),N=norm(M.cross(T,U));
  if(L<1e-12||width<1e-12)return {point:base,invalidFrame:true,stretch:null,determinant:null,reason:'singular-posed-chart'};
  const sc=root.VolumeConstraint.scales(L,width,B.sheet,root.TendonDeformation.fraction(w,B.definition.tissues)),off=sub(p,r.c),candidate=add(s.c,add(add(mul(T,dot(off,R.t)*sc.t),mul(add(mul(T,shear),mul(U,sc.u)),dot(off,R.u))),mul(N,dot(off,R.n)*sc.n))),gate=M.smooth(w/.12)*M.smooth((1-w)/.12),point=root.ContactConstraint.project(p,mix(base,candidate,gate),w,state.contact);
  return {point,stretch:sc.t,determinant:sc.determinant,bounded:sc.bounded,chartSpacingMm:sigma/B.definition.source.worldUnitsPerMm};
 }
 function pose(B,a,b,contact=null,originAt=null){const s=E.pose(B,a,b,contact,originAt);s.vertexDiagnostics=i=>evaluate(s,i);return s;}
 function deform(state,out=new Float32Array(state.binding.rest.length)){for(let i=0;i<out.length/3;i++)out.set(evaluate(state,i).point,i*3);return out;}
 root.DifferentialChart=Object.freeze({prepare,pose,evaluate,deform});
})(globalThis);
