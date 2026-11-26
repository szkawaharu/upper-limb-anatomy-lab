/* Deterministic correction field in joint coordinates.
 * No temporal lag, muscle-specific offsets, or change to joint commands.
 * Interpolation is not proof of contact feasibility: audit at query pose.
 */
(function(root){
"use strict";
const M=root.AnatomyMath, I=[0,0,0,1];
const axes=[[0,180,15],[0,90,15],[-90,90,45]];
function canon(q){const n=Math.hypot(...q);if(!Number.isFinite(n)||n<1e-12)throw Error("invalid rotation");return q.map(x=>x*(q[3]<0?-1:1)/n)}
function correction(raw,solved){return canon(M.qmul(solved,[-raw[0],-raw[1],-raw[2],raw[3]]))}
function bracket(x,[lo,hi,step]){
 if(!Number.isFinite(x)||x<lo-1e-8||x>hi+1e-8)throw Error("outside pose-field domain");
 x=M.clamp(x,lo,hi);const a=Math.min(hi,lo+Math.floor((x-lo)/step)*step),b=Math.min(hi,a+step);
 if(a===b||x===a)return [[a,1]];
 if(x===b)return [[b,1]];
 const t=(x-a)/(b-a);return [[a,1-t],[b,t]];
}
function create(solveNode,grid=axes){
 const cache=new Map();let coldMs=0,coldMaxMs=0;
 function node(e,p,a){
  const k=[e,p,a].join(":");
  if(!cache.has(k)){
   const t=performance.now(),n=solveNode(e,p,a);
   if(!n?.raw||!n?.solved)throw Error("pose node unavailable");
   cache.set(k,{C:correction(n.raw.qC,n.solved.qC),S:correction(n.raw.qS,n.solved.qS)});
   const dt=performance.now()-t;coldMs+=dt;coldMaxMs=Math.max(coldMaxMs,dt);
  }
  return cache.get(k);
 }
 function sample(e,p,a,raw,SC,cv){
  const C=[0,0,0,0],S=[0,0,0,0];
  for(const [x,u] of bracket(e,grid[0]))for(const [y,v] of bracket(p,grid[1]))for(const [z,w] of bracket(a,grid[2])){
   const n=node(x,y,z),f=u*v*w;
   for(let i=0;i<4;i++){C[i]+=f*n.C[i];S[i]+=f*n.S[i]}
  }
  const qC=canon(M.qmul(canon(C),raw.qC)),qS=canon(M.qmul(canon(S),raw.qS));
  return {qC,qS,AC:M.add(SC,M.qrot(qC,cv))};
 }
 return {sample,stats:()=>({nodes:cache.size,coldMs,coldMaxMs}),grid};
}
root.ShoulderPoseField=Object.freeze({create,correction,axes});
})(globalThis);
