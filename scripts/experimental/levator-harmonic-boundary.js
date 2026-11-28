/* Experimental five-bone harmonic coordinate control, not physical elasticity.
 * Positive inverse-length graph energy, exact Dirichlet values. Components
 * without any pins retain the explicitly reported legacy two-bone base map.
 * Never remove source triangles or impose fictitious attachment points. */
(function(root){
'use strict';
const M=root.AnatomyMath;
function prepare(data){
 const reg=root.LevatorBoundaryRegistry.registry(data),n=reg.label.length,adj=Array.from({length:n},()=>new Map()),rest=data.rest;
 for(let k=0;k<data.faces.length;k+=3)for(let j=0;j<3;j++){
  const a=data.faces[k+j],b=data.faces[k+(j+1)%3];
  if(a===b||a<0||b<0||a>=n||b>=n)throw Error('Invalid source edge');
  const d=Math.hypot(...[0,1,2].map(i=>rest[3*a+i]-rest[3*b+i]));
  if(!(d>0))throw Error('Zero source edge');
  adj[a].set(b,1/d);adj[b].set(a,1/d);
 }
 const seen=new Uint8Array(n),weights=new Float64Array(n*5),unanchored=new Uint8Array(n),components=[];
 for(let i=0;i<n;i++)if(!seen[i]){
  const ids=[i];seen[i]=1;
  for(let j=0;j<ids.length;j++)for(const k of adj[ids[j]].keys())if(!seen[k]){seen[k]=1;ids.push(k)}
  const pins=ids.filter(j=>reg.label[j]>=0),labels=[...new Set(pins.map(j=>reg.label[j]))],free=ids.filter(j=>reg.label[j]<0);
  const info={vertices:ids.length,pins:pins.length,labels:labels.map(j=>reg.keys[j]),iterations:[],maxRelativeResidual:0};
  components.push(info);
  if(!labels.length){for(const j of ids)unanchored[j]=1;continue}
  if(labels.length===1){for(const j of ids)weights[5*j+labels[0]]=1;continue}
  for(const j of pins)weights[5*j+reg.label[j]]=1;
  if(!free.length)continue;
  const index=new Map(free.map((j,k)=>[j,k])),diag=Float64Array.from(free,j=>[...adj[j].values()].reduce((a,b)=>a+b,0));
  const A=x=>Float64Array.from(free,(j,k)=>{let v=diag[k]*x[k];for(const [l,w]of adj[j])if(index.has(l))v-=w*x[index.get(l)];return v});
  const dot=(x,y)=>x.reduce((s,v,j)=>s+v*y[j],0);
  for(const label of labels){
   const rhs=Float64Array.from(free,j=>[...adj[j]].reduce((s,[k,w])=>s+(reg.label[k]===label?w:0),0));
   const x=new Float64Array(free.length),r=rhs.slice(),z=Float64Array.from(r,(v,j)=>v/diag[j]),p=z.slice(),scale=Math.max(Math.sqrt(dot(rhs,rhs)),1e-30);
   let rz=dot(r,z),it=0,res=Math.sqrt(dot(r,r))/scale;
   for(;res>1e-10&&it<Math.max(1000,free.length*2);it++){
    const Ap=A(p),den=dot(p,Ap);
    if(!(den>0))throw Error('Nonpositive graph metric');
    const a=rz/den;
    for(let j=0;j<x.length;j++){x[j]+=a*p[j];r[j]-=a*Ap[j];z[j]=r[j]/diag[j]}
    const next=dot(r,z),beta=next/rz;
    for(let j=0;j<p.length;j++)p[j]=z[j]+beta*p[j];
    rz=next;res=Math.sqrt(dot(r,r))/scale;
   }
   const actual=A(x);res=Math.sqrt(actual.reduce((s,v,j)=>s+(v-rhs[j])**2,0))/scale;
   if(res>1e-9)throw Error('Harmonic solve failed residual '+res);
   info.iterations.push(it);info.maxRelativeResidual=Math.max(info.maxRelativeResidual,res);
   for(let k=0;k<free.length;k++)weights[5*free[k]+label]=x[k];
  }
 }
 let partitionError=0,min=1,max=0;
 for(let i=0;i<n;i++)if(!unanchored[i]){
  const a=weights.slice(i*5,i*5+5),sum=a.reduce((s,v)=>s+v,0);
  partitionError=Math.max(partitionError,Math.abs(sum-1));min=Math.min(min,...a);max=Math.max(max,...a);
  if(Math.abs(sum-1)>1e-7||Math.min(...a)<-1e-7||Math.max(...a)>1+1e-7)throw Error('Invalid harmonic partition');
  for(let j=0;j<5;j++)weights[i*5+j]=Math.max(0,a[j])/a.reduce((s,v)=>s+Math.max(0,v),0);
 }
 return {reg,weights,unanchored,diagnostics:{components,unanchoredVertices:unanchored.reduce((a,b)=>a+b,0),partitionError,min,max,physicalVolumeVerified:false}};
}
function deform(state,prepared,bones){
 const B=state.binding,{reg,weights,unanchored}=prepared,out=new Float32Array(B.rest.length);
 for(let i=0;i<B.weights.length;i++){
  const p=Array.from(B.rest.slice(i*3,i*3+3));let x;
  if(unanchored[i]){
   const a=state.originAt?state.originAt(p):state.a;x=M.mix(M.apply(a,p),M.apply(state.b,p),B.weights[i]);
  }else if(reg.label[i]>=0)x=M.apply(bones[reg.keys[reg.label[i]]],p);
  else {
   x=[0,0,0];
   for(let j=0;j<5;j++)if(weights[i*5+j])x=M.add(x,M.mul(M.apply(bones[reg.keys[j]],p),weights[i*5+j]));
  }
  if(!x.every(Number.isFinite))throw Error('Nonfinite surface target');
  out.set(x,3*i);
 }
 return out;
}
root.LevatorHarmonicBoundary=Object.freeze({prepare,deform});
})(globalThis);
