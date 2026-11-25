/* ScapulothoracicEngineV066
 * Rigid scapula constrained against a thoracic surface proxy in three regions.
 * Geometry-first educational/research model. No cartilage, bursa, contact force,
 * ligament, muscle-force, pathology or clinical-normality solver.
 */
(function(root){
'use strict';
const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
const smooth=x=>{x=clamp(x,0,1);return x*x*(3-2*x)};
const add=(a,b)=>a.map((x,i)=>x+b[i]),sub=(a,b)=>a.map((x,i)=>x-b[i]),mul=(a,s)=>a.map(x=>x*s);
const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0),cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const len=a=>Math.hypot(...a),norm=a=>mul(a,1/Math.max(len(a),1e-12));
const qa=(axis,a)=>{const n=norm(axis),s=Math.sin(a/2);return [n[0]*s,n[1]*s,n[2]*s,Math.cos(a/2)]};
const qm=(a,b)=>[a[3]*b[0]+a[0]*b[3]+a[1]*b[2]-a[2]*b[1],a[3]*b[1]-a[0]*b[2]+a[1]*b[3]+a[2]*b[0],a[3]*b[2]+a[0]*b[1]-a[1]*b[0]+a[2]*b[3],a[3]*b[3]-a[0]*b[0]-a[1]*b[1]-a[2]*b[2]];
const qrot=(q,v)=>{const u=q.slice(0,3),t=mul(cross(u,v),2);return add(v,add(mul(t,q[3]),cross(u,t)))};
const transform=(q,o,p)=>add(o,qrot(q,p));
const rad=d=>d*Math.PI/180;

function regionOf(p,lo,hi){
 const y=clamp((p[1]-lo[1])/Math.max(hi[1]-lo[1],1e-9),0,1);
 return y>.66?'upper':y>.34?'middle':'lower';
}
function build(points,restAC,profileNormal){
 if(!points?.length||!restAC||typeof profileNormal!=='function')return null;
 const lo=[Infinity,Infinity,Infinity],hi=[-Infinity,-Infinity,-Infinity];
 for(const v of points)for(let i=0;i<3;i++){lo[i]=Math.min(lo[i],v.p[i]);hi[i]=Math.max(hi[i],v.p[i])}
 const groups={upper:[],middle:[],lower:[]};
 for(const v of points){const region=v.region||regionOf(v.p,lo,hi);groups[region].push({...v,region})}
 if(Object.values(groups).some(x=>x.length<2))return null;
 const centroid={};
 for(const [name,arr] of Object.entries(groups)){
  let c=[0,0,0],w=0;for(const v of arr){const k=v.weight||1;c=add(c,mul(add(restAC,v.p),k));w+=k}centroid[name]=mul(c,1/Math.max(w,1e-9));
 }
 // A stable whole-blade costal normal: medial-border superior->inferior crossed
 // with medial->lateral direction. The sign is irrelevant; tangency uses abs(dot).
 const upper=centroid.upper,lower=centroid.lower,mid=centroid.middle;
 const medial=norm(sub(lower,upper));
 const lateralHint=norm([1,0,0]);
 let restNormal=norm(cross(medial,lateralHint));
 if(len(restNormal)<.5)restNormal=[0,0,1];
 const restTangency={};
 for(const name of ['upper','middle','lower']){
  const n=profileNormal(centroid[name]);restTangency[name]=n?Math.abs(dot(restNormal,n)):null;
 }
 return {groups,centroid,restNormal,restTangency,lo,hi};
}
function metrics(model,qS,AC,gapFn,normalFn,scale){
 const regions={},all=[];let globalW=0,globalSq=0,maxOut=-Infinity,maxIn=Infinity;
 for(const name of ['upper','middle','lower']){
  const arr=model.groups[name];let w=0,sq=0,out=-Infinity,inside=Infinity,n=0;
  for(const v of arr){
   if(v.sliceRest==null)continue;
   const p=transform(qS,AC,v.p),g=gapFn(p);if(g==null)continue;
   const d=(g-v.sliceRest)/scale,k=v.weight||1;w+=k;sq+=k*d*d;out=Math.max(out,d);inside=Math.min(inside,d);n++;all.push({d,k});
  }
  const c=transform(qS,AC,sub(model.centroid[name],AC)); // only used if AC is rest-like; replaced below by local centroid calc
  // exact transformed centroid using rest local coordinates:
  const local=sub(model.centroid[name],model.restAC||AC),world=transform(qS,AC,local),tn=normalFn(world);
  const sn=qrot(qS,model.restNormal),tangency=tn?Math.abs(dot(norm(sn),norm(tn))):null,rest=model.restTangency[name];
  regions[name]={n,rms:w?Math.sqrt(sq/w):Infinity,maxOut:out,maxIn:inside,tangency,tangencyDrop:tangency!=null&&rest!=null?Math.max(0,rest-tangency):null};
  globalW+=w;globalSq+=sq;maxOut=Math.max(maxOut,out);maxIn=Math.min(maxIn,inside);
 }
 return {regions,rms:globalW?Math.sqrt(globalSq/globalW):Infinity,maxOut,maxIn};
}
function score(m){
 if(!Number.isFinite(m.rms))return Infinity;
 const cfg={
  upper:{w:1.0,out:10,inward:5,tan:1.0},
  middle:{w:1.35,out:8,inward:4,tan:1.4},
  lower:{w:2.15,out:6,inward:4,tan:2.1}
 };
 let c=.55*m.rms*m.rms;
 for(const name of ['upper','middle','lower']){
  const r=m.regions[name],g=cfg[name];if(!r||!Number.isFinite(r.rms))return Infinity;
  c+=g.w*r.rms*r.rms;
  c+=g.w*3.2*Math.max(0,r.maxOut-g.out)**2;
  c+=g.w*4.3*Math.max(0,-r.maxIn-g.inward)**2;
  if(Number.isFinite(r.tangencyDrop))c+=g.tan*900*r.tangencyDrop*r.tangencyDrop;
 }
 return c;
}
// Deterministic bounded coordinate line minimization. No pose-history state.
// Preserve the supplied objective and hard feasibility checks; never adopt a
// non-finite or uphill trial. Continuous-valued trial angles replace a fixed lattice.
// This is an experimental local optimizer, not a proof of global continuity.
function minimize(evaluate,limits,passes=4){
 let a=limits.map(()=>0),best=evaluate(a);
 for(let sweep=0;sweep<passes;sweep++){
  for(let axis=0;axis<a.length;axis++){
   const L=limits[axis];if(L<1e-10)continue;
   const from=a.slice();
   const at=x=>{const p=from.slice();p[axis]=x;return {p,r:evaluate(p)}};
   const accept=o=>{if(o.r.total<best.total-1e-8){a=o.p;best=o.r}};
   // Survey the complete bounded line so the bracket is not tied to the
   // previous pose or the sign/order of greedy trial steps.
   const xs=Array.from({length:9},(_,i)=>-L+2*L*i/8).concat(from[axis]).sort((x,y)=>x-y);
   const vals=xs.map(at);let j=0;for(let i=1;i<vals.length;i++)if(vals[i].r.total<vals[j].r.total)j=i;
   for(const o of vals)accept(o);
   let lo=xs[Math.max(0,j-1)],hi=xs[Math.min(xs.length-1,j+1)];
   // Golden-section line refinement, including feasible incumbent.
   for(let k=0;k<10&&hi-lo>1e-6;k++){
    const x=lo+(hi-lo)*.3819660112501051,y=hi-(hi-lo)*.3819660112501051;
    const u=at(x),v=at(y);accept(u);accept(v);
    if(u.r.total<=v.r.total)hi=y;else lo=x;
   }
  }
 }
 return {a,best};
}
function solve(opts){
 const {baseQC,baseQS,SC,cv,restAC,elev=0,plane=0,points,gapFn,normalFn,scale=1}=opts;
 const identity={qC:baseQC,qS:baseQS,AC:add(SC,qrot(baseQC,cv)),before:null,after:null,params:[0,0,0,0,0,0],metrics:null};
 if(!points?.length||typeof gapFn!=='function'||typeof normalFn!=='function')return identity;
 const model=build(points,restAC,normalFn);if(!model)return identity;model.restAC=restAC.slice();
 const baseAC=add(SC,qrot(baseQC,cv)),baseM=metrics(model,baseQS,baseAC,gapFn,normalFn,scale),baseScore=score(baseM);
 if(!Number.isFinite(baseScore))return {...identity,metrics:baseM};
 const k=smooth(clamp(elev/45,0,1)),limits=[24*k,24*k,20*k,7*k,9*k,8*k];
 const make=a=>{
  const qC0=qm(qa([0,0,1],rad(a[3])),qm(qa([0,1,0],rad(a[4])),baseQC));
  const axis=norm(qrot(qC0,cv)),qC=qm(qa(axis,rad(a[5])),qC0),AC=add(SC,qrot(qC,cv));
  const qS=qm(qa([0,0,1],rad(a[2])),qm(qa([0,1,0],rad(a[0])),qm(qa([1,0,0],rad(a[1])),baseQS)));
  return {qC,qS,AC};
 };
 const evaluate=a=>{
  const p=make(a),m=metrics(model,p.qS,p.AC,gapFn,normalFn,scale),raw=score(m);
  // Corrections remain small and continuous; geometry dominates. Plane-dependent
  // regularizer is deliberately weak because normal kinematics vary substantially.
  const reg=.10*(a[0]**2+a[1]**2)+.08*a[2]**2+.75*(a[3]**2+a[4]**2)+.55*a[5]**2+.002*Math.abs(plane)*a[0]*a[0];
  return {p,m,raw,total:raw+reg};
 };
 const {a,best}=minimize(evaluate,limits,4);
 // Never adopt a solution that worsens the geometry-only objective.
 if(best.raw>baseScore+1e-8)return {...identity,before:baseScore,after:baseScore,metrics:baseM};
 return {...best.p,before:baseScore,after:best.raw,params:a,metrics:best.m,model};
}
function humeralAudit(qS,qH,restShaft){
 const shaft=norm(qrot(qH,restShaft)),local=norm(qrot([-qS[0],-qS[1],-qS[2],qS[3]],shaft));
 // Descriptive only: angle to scapular-plane normal proxy.
 return {shaftWorld:shaft,shaftInScapula:local,planeNormalAngleDeg:Math.acos(clamp(Math.abs(local[2]),-1,1))*180/Math.PI};
}
root.ScapulothoracicEngineV066=Object.freeze({version:'0.66.0',build,metrics,score,solve,humeralAudit,regionOf,minimize});
})(globalThis);
