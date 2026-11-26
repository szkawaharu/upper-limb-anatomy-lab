/* Conservative source-bone support enclosure. For halfspaces d_i <= 0,
 * tau*log(mean(exp(d_i/tau))) <= max(d_i). Thus a nonnegative gap implies
 * separation from the enclosing polytope, and therefore from its source ribs.
 * This support value is NOT Euclidean distance or measured fascia thickness.
 */
(function(root){
 'use strict';
 const M=root.AnatomyMath;
 function fromAtlas(data,SC,atlasSC,scale){
  const map=p=>[-p[0],p[2],p[1]],offset=M.sub(M.mul(SC,1/scale),map(atlasSC));
  return {normals:data.normals.map(map),supportsMm:data.supportsMm.map((d,i)=>d+M.dot(map(data.normals[i]),offset)),points:data.scapularPointsMm.map(p=>M.mul(M.add(map(p),offset),scale)),smoothingMm:data.smoothingMm,residualScaleMm:data.residualScaleMm,sourceCommit:data.sourceCommit,sourceMethod:data.method};
 }
 const FREE=Object.freeze({gap:0,gradient:Object.freeze([0,0,0])});
 function evaluate(p,body,scale,gradient=false,negativeOnly=false){
  const x=p[0]/scale,y=p[1]/scale,z=p[2]/scale,ns=body.normals,ds=body.supportsMm,tau=body.smoothingMm;
  let max=-Infinity;
  const safe=tau*Math.log(ns.length);
  // log-mean-exp >= max(d)-tau*log(N). One sufficiently separating
  // halfspace therefore proves a zero unilateral residual without evaluating
  // every exponential. This is an exact broad phase, not a distance cutoff.
  for(let i=0;i<ns.length;i++){const n=ns[i],d=n[0]*x+n[1]*y+n[2]*z-ds[i];if(negativeOnly&&d>=safe)return gradient?FREE:0;if(d>max)max=d;}
  let sum=0,gx=0,gy=0,gz=0;
  for(let i=0;i<ns.length;i++){const n=ns[i],w=Math.exp((n[0]*x+n[1]*y+n[2]*z-ds[i]-max)/tau);sum+=w;if(gradient){gx+=w*n[0];gy+=w*n[1];gz+=w*n[2];}}
  const gap=max+tau*Math.log(sum/ns.length);
  return gradient?{gap,gradient:[gx/sum,gy/sum,gz/sum]}:gap;
 }
 const gapMm=(p,body,scale)=>evaluate(p,body,scale,false);
 const penetration=(p,body,scale,gradient=false)=>evaluate(p,body,scale,gradient,true);
 root.ConvexContact=Object.freeze({fromAtlas,gapMm,evaluate,penetration});
})(globalThis);
