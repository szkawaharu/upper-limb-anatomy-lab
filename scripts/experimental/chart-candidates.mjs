// Compare shared chart constructions on a captured, immutable atlas binding.
// No visibility, angle, contact, or per-muscle coefficient changes.
export function rebuildChart(B,kind) {
 const M=globalThis.AnatomyMath;
 if(kind==='existing')return B;
 const controls=B.field.controls.map(c=>{
  if(kind==='linear')return [c[0],M.mix(c[0],c[3],1/3),M.mix(c[0],c[3],2/3),c[3]];
  return c;
 });
 return {...B,field:{...B.field,controls},material:kind==='linear'?B.weights:B.material};
}
export function restSurfacePoint(state,p,w,l,t,kind){
 const M=AnatomyMath,P=FiberPath,{mix,apply,sub,add,mul,qrot,qmul,qblendLift,dot,norm,len}=M,{binding:B,a,b}=state;
 if(w<=0||w>=1)return apply(w>=1?b:a,p);
 const q=qblendLift(a.q,b.q,w);
 if(kind==='dual-quat'){
  const dual=A=>qmul([...sub(A.p,qrot(A.q,A.rest)),0],A.q).map(x=>x*.5),da=dual(a),db=dual(b),raw=a.q.map((v,i)=>v*(1-w)+b.q[i]*w),scale=Math.hypot(...raw),d=da.map((v,i)=>(v*(1-w)+db[i]*w)/scale),v=qmul(d,[-q[0],-q[1],-q[2],q[3]]);
  return add(qrot(q,p),v.slice(0,3).map(x=>2*x));
 }
 const c=P.laneCurve(B.field.controls,l),C=P.curve(c,t),off=sub(p,C),center=mix(apply(a,C),apply(b,C),w);
 if(kind==='center')return add(center,qrot(q,off));
 const T=norm(sub(c[3],c[0])),ref=sub(B.reference,mul(T,dot(T,B.reference))),U=norm(ref),N=norm(M.cross(T,U)),lambda=len(sub(apply(b,c[3]),apply(a,c[0])))/len(sub(c[3],c[0])),tendon=TendonDeformation.fraction(w,B.definition.tissues),sc=VolumeConstraint.scales(lambda,1,B.sheet,tendon),scaled=add(mul(T,dot(off,T)),add(mul(U,dot(off,U)*sc.u),mul(N,dot(off,N)*sc.n)));
 return add(center,qrot(q,scaled));
}
