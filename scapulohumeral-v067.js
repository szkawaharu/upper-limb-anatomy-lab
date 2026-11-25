/* ScapulohumeralGuardV067
 * Literature-informed broad kinematic guard for the v0.66 scapulothoracic solver.
 * Proxy angles are atlas-coordinate Y-Z-X Euler components, NOT clinical ISB angles.
 * The guard preserves commanded humeral orientation and only limits excessive
 * scapular/clavicular correction when thorax geometry and GH contribution remain safe.
 */
(function(root){
'use strict';
const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
const add=(a,b)=>a.map((x,i)=>x+b[i]),sub=(a,b)=>a.map((x,i)=>x-b[i]),mul=(a,s)=>a.map(x=>x*s);
const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0),cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const len=a=>Math.hypot(...a),norm=a=>mul(a,1/Math.max(len(a),1e-12));
const qnorm=q=>{const n=Math.hypot(...q)||1;return q.map(x=>x/n)};
const qconj=q=>[-q[0],-q[1],-q[2],q[3]];
const qmul=(a,b)=>[a[3]*b[0]+a[0]*b[3]+a[1]*b[2]-a[2]*b[1],a[3]*b[1]-a[0]*b[2]+a[1]*b[3]+a[2]*b[0],a[3]*b[2]+a[0]*b[1]-a[1]*b[0]+a[2]*b[3],a[3]*b[3]-a[0]*b[0]-a[1]*b[1]-a[2]*b[2]];
const qrot=(q,v)=>{const u=q.slice(0,3),t=mul(cross(u,v),2);return add(v,add(mul(t,q[3]),cross(u,t)))};
const deg=r=>r*180/Math.PI;
function slerp(a,b,t){
 a=qnorm(a);b=qnorm(b);let d=dot(a,b);if(d<0){b=mul(b,-1);d=-d}
 if(d>.9995)return qnorm(a.map((x,i)=>x+(b[i]-x)*t));
 const th=Math.acos(clamp(d,-1,1)),s=Math.sin(th);return a.map((x,i)=>(Math.sin((1-t)*th)*x+Math.sin(t*th)*b[i])/s);
}
function qangle(a,b){let d=Math.abs(dot(qnorm(a),qnorm(b)));return deg(2*Math.acos(clamp(d,-1,1)))}
function matrix(q){
 q=qnorm(q);const [x,y,z,w]=q,xx=x*x,yy=y*y,zz=z*z,xy=x*y,xz=x*z,yz=y*z,wx=w*x,wy=w*y,wz=w*z;
 return [[1-2*(yy+zz),2*(xy-wz),2*(xz+wy)],[2*(xy+wz),1-2*(xx+zz),2*(yz-wx)],[2*(xz-wy),2*(yz+wx),1-2*(xx+yy)]];
}
// Decompose the atlas proxy as Ry(y) * Rz(z) * Rx(x).
function scapProxy(q){
 const R=matrix(q),z=Math.asin(clamp(R[1][0],-1,1)),cz=Math.cos(z);
 if(Math.abs(cz)<1e-5)return {internalDeg:deg(Math.atan2(R[0][2],R[2][2])),upwardDeg:deg(z),tiltDeg:0,singular:true};
 return {internalDeg:deg(Math.atan2(-R[2][0],R[0][0])),upwardDeg:deg(z),tiltDeg:deg(Math.atan2(-R[1][2],R[1][1])),singular:false};
}
function envelope(elev,plane){
 const k=clamp(elev/180,0,1),forward=clamp(Math.sin(plane*Math.PI/180),-1,1);
 const center={upwardDeg:elev/3,tiltDeg:13*k,internalDeg:(-8+12*forward)*k};
 return {
  center,
  tol:{upwardDeg:10+5*k,tiltDeg:8+5*k,internalDeg:14+4*k},
  maxCorrectionDeg:10+10*k,
  // Broad descriptive GH contribution band. It is intentionally much wider
  // than a single scapulohumeral-rhythm ratio.
  ghMinDeg:Math.max(0,elev*(.46+.02*forward)-15),
  ghMaxDeg:elev*(.90+.02*forward)+15,
  ghDeviationFromBaseMaxDeg:12+8*k
 };
}
function relativeGH(qS,qH,restShaft){
 const rel=qmul(qconj(qnorm(qS)),qnorm(qH)),axis=norm(restShaft),moved=norm(qrot(rel,axis));
 return deg(Math.acos(clamp(dot(axis,moved),-1,1)));
}
function withinProxy(p,e){
 return Math.abs(p.upwardDeg-e.center.upwardDeg)<=e.tol.upwardDeg &&
        Math.abs(p.tiltDeg-e.center.tiltDeg)<=e.tol.tiltDeg &&
        Math.abs(p.internalDeg-e.center.internalDeg)<=e.tol.internalDeg;
}
function audit(qS,qH,restShaft,elev,plane,baseQS=qS){
 const env=envelope(elev,plane),proxy=scapProxy(qS),baseProxy=scapProxy(baseQS),gh=relativeGH(qS,qH,restShaft),baseGH=relativeGH(baseQS,qH,restShaft),correctionDeg=qangle(baseQS,qS);
 return {proxy,baseProxy,env,ghDeg:gh,baseGHDeg:baseGH,correctionDeg,
  proxyOK:withinProxy(proxy,env),
  ghOK:gh>=env.ghMinDeg&&gh<=env.ghMaxDeg&&Math.abs(gh-baseGH)<=env.ghDeviationFromBaseMaxDeg,
  correctionOK:correctionDeg<=env.maxCorrectionDeg};
}
function refine(opts){
 const {v66,base,candidate,SC,cv,restAC,elev=0,plane=0,qH,restShaft,points,gapFn,normalFn,scale=1}=opts;
 if(!v66||!candidate||!qH||!restShaft)return {...candidate,fallback:'v0.66',v67:null};
 const model=v66.build(points,restAC,normalFn);if(!model)return {...candidate,fallback:'v0.66',v67:null};
 model.restAC=restAC.slice();
 const geom=(qS,AC)=>v66.score(v66.metrics(model,qS,AC,gapFn,normalFn,scale));
 const baseGeom=geom(base.qS,base.AC),candidateGeom=geom(candidate.qS,candidate.AC),baseAudit=audit(base.qS,qH,restShaft,elev,plane,base.qS);
 let chosen={qC:base.qC,qS:base.qS,AC:base.AC,t:0,geom:baseGeom,audit:baseAudit};
 // Keep as much of the v0.66 correction as possible, but do not accept
 // a kinematic guard violation or geometry worse than the uncorrected pose.
 for(const t of [1,.875,.75,.625,.5,.375,.25,.125]){
  const qC=slerp(base.qC,candidate.qC,t),qS=slerp(base.qS,candidate.qS,t),AC=add(SC,qrot(qC,cv)),g=geom(qS,AC),a=audit(qS,qH,restShaft,elev,plane,base.qS);
  if(g<=baseGeom+1e-8&&a.proxyOK&&a.ghOK&&a.correctionOK){chosen={qC,qS,AC,t,geom:g,audit:a};break}
 }
 return {...chosen,before:baseGeom,after:chosen.geom,v66Before:candidate.before,v66After:candidate.after,
  v66CandidateGeom:candidateGeom,v67:{blend:chosen.t,...chosen.audit},fallback:chosen.t===0?'base-guard':null};
}
root.ScapulohumeralGuardV067=Object.freeze({version:'0.67.0',scapProxy,envelope,relativeGH,audit,refine,slerp,qangle});
})(globalThis);
