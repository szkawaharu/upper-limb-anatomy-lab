/* ScapulothoracicContactGuardV072
 * Contact-priority post-refinement for the scapula on the thorax.
 *
 * v0.66 optimizes broad costal-surface fit; v0.67 can deliberately blend that
 * correction back toward the literature-informed base pose. v0.72 runs AFTER
 * v0.67 and intervenes only when local costal-surface separation/penetration
 * exceeds engineering guard thresholds. It keeps the commanded humeral pose
 * unchanged and searches only small scapular/clavicular corrections.
 *
 * Thresholds are numerical safeguards, NOT clinical normal/pathologic limits.
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
const rad=d=>d*Math.PI/180;

const CFG=Object.freeze({
 upper:{w:1.0,out:10,inward:5,rms:8,tan:.22},
 middle:{w:1.55,out:8,inward:4.5,rms:6.5,tan:.18},
 lower:{w:2.8,out:6,inward:4.5,rms:5.5,tan:.14}
});

function regionPenalty(r,c){
 if(!r||!Number.isFinite(r.rms))return Infinity;
 const out=Math.max(0,r.maxOut-c.out),inside=Math.max(0,-r.maxIn-c.inward),rms=Math.max(0,r.rms-c.rms);
 const tan=Number.isFinite(r.tangencyDrop)?Math.max(0,r.tangencyDrop-c.tan):0;
 return c.w*(10000*(r.missing||0)/Math.max(1,r.expected||r.n||1)+.12*r.rms*r.rms+7.5*out*out+10*inside*inside+2.0*rms*rms+1500*tan*tan);
}
function violation(r,c){
 if(!r||!Number.isFinite(r.rms))return {out:Infinity,inside:Infinity,rms:Infinity,tan:Infinity};
 return {
  out:Math.max(0,r.maxOut-c.out),
  inside:Math.max(0,-r.maxIn-c.inward),
  rms:Math.max(0,r.rms-c.rms),
  tan:Number.isFinite(r.tangencyDrop)?Math.max(0,r.tangencyDrop-c.tan):0
 };
}
function noRegionalRegression(m,base){
 if(!m||!base)return false;
 for(const k of ['upper','middle','lower']){
  if((m.regions[k]?.missing||0)>(base.regions[k]?.missing||0))return false;
  const a=violation(m.regions[k],CFG[k]),b=violation(base.regions[k],CFG[k]);
  if(a.out>b.out+.5||a.inside>b.inside+.5||a.rms>b.rms+.5||a.tan>b.tan+.02)return false;
 }
 return true;
}
function contactScore(m){
 if(!m||!Number.isFinite(m.rms))return Infinity;
 let s=.12*m.rms*m.rms;
 for(const k of ['upper','middle','lower'])s+=regionPenalty(m.regions[k],CFG[k]);
 const globalOut=Math.max(0,m.maxOut-10),globalIn=Math.max(0,-m.maxIn-5);
 return s+5*globalOut*globalOut+7*globalIn*globalIn;
}
function contactOK(m){
 if(!m||!Number.isFinite(m.rms))return false;
 for(const k of ['upper','middle','lower']){
  const r=m.regions[k],c=CFG[k];if(!r||!Number.isFinite(r.rms))return false;
  if((r.missing||0)>0||r.maxOut>c.out||r.maxIn<-c.inward||r.rms>c.rms)return false;
  if(Number.isFinite(r.tangencyDrop)&&r.tangencyDrop>c.tan)return false;
 }
 return true;
}
function maxOut(m){return Math.max(...['upper','middle','lower'].map(k=>m?.regions?.[k]?.maxOut??-Infinity))}
function lowerOut(m){return m?.regions?.lower?.maxOut??Infinity}

function hardKinematicOK(a){
 if(!a?.env||!a?.proxy||!Number.isFinite(a.correctionDeg))return true;
 const e=a.env,p=a.proxy;
 const proxy=
  Math.abs(p.upwardDeg-e.center.upwardDeg)<=e.tol.upwardDeg+8 &&
  Math.abs(p.tiltDeg-e.center.tiltDeg)<=e.tol.tiltDeg+8 &&
  Math.abs(p.internalDeg-e.center.internalDeg)<=e.tol.internalDeg+10;
 const gh=a.ghDeg>=e.ghMinDeg-10 && a.ghDeg<=e.ghMaxDeg+10 &&
  Math.abs(a.ghDeg-a.baseGHDeg)<=e.ghDeviationFromBaseMaxDeg+10;
 return proxy&&gh&&a.correctionDeg<=e.maxCorrectionDeg+8;
}

function refine(opts){
 const {v66,v67,input,base,SC,cv,restAC,elev=0,plane=0,qH,restShaft,points,gapFn,normalFn,scale=1}=opts||{};
 if(!v66||!input||!base||!SC||!cv||!restAC||!points?.length||typeof gapFn!=='function'||typeof normalFn!=='function'){
  return {...input,fallback072:'missing-input',v72:null};
 }
 const model=v66.build(points,restAC,normalFn);
 if(!model)return {...input,fallback072:'model-unavailable',v72:null};
 model.restAC=restAC.slice();
 const metrics=(qS,AC)=>v66.metrics(model,qS,AC,gapFn,normalFn,scale);
 const beforeM=metrics(input.qS,input.AC),beforeScore=contactScore(beforeM);
 if(!Number.isFinite(beforeScore))return {...input,fallback072:'metrics-unavailable',v72:null};
 if(contactOK(beforeM)){
  return {...input,v72:{active:false,reason:'contact-within-guard',before:beforeM,after:beforeM,beforeScore,afterScore:beforeScore,params:[0,0,0,0,0,0],correctionDeg:0}};
 }

 const k=smooth(clamp(elev/45,0,1));
 // Small post-polish only. The broad motion remains controlled by v0.66/v0.67.
 const limits=[12*k,14*k,10*k,6*k,8*k,6*k];
 const make=a=>{
  const qC0=qm(qa([0,0,1],rad(a[3])),qm(qa([0,1,0],rad(a[4])),input.qC));
  const axis=norm(qrot(qC0,cv)),qC=qm(qa(axis,rad(a[5])),qC0),AC=add(SC,qrot(qC,cv));
  const qS=qm(qa([0,0,1],rad(a[2])),qm(qa([0,1,0],rad(a[0])),qm(qa([1,0,0],rad(a[1])),input.qS)));
  return {qC,qS,AC};
 };
 const evaluate=a=>{
  const p=make(a),m=metrics(p.qS,p.AC),raw=contactScore(m);
  const audit=v67&&qH&&restShaft?v67.audit(p.qS,qH,restShaft,elev,plane,base.qS):null;
  if(audit&&!hardKinematicOK(audit))return {p,m,raw,total:Infinity,audit};
  if(!noRegionalRegression(m,beforeM))return {p,m,raw,total:Infinity,audit};
  const reg=.20*(a[0]**2+a[1]**2)+.16*a[2]**2+.85*(a[3]**2+a[4]**2)+.65*a[5]**2;
  return {p,m,raw,total:raw+reg,audit};
 };

 let a=[0,0,0,0,0,0],best=evaluate(a);
 for(const step of [4,2,1,.5,.25]){
  let changed=true,loops=0;
  while(changed&&loops++<3){
   changed=false;
   for(let axis=0;axis<6;axis++){
    const from=a.slice();
    for(const sign of [-1,1]){
     const c=from.slice();c[axis]=clamp(c[axis]+sign*step,-limits[axis],limits[axis]);
     const t=evaluate(c);
     if(t.total<best.total-1e-7){a=c;best=t;changed=true}
    }
   }
  }
 }

 // v0.75: coordinated scapula/clavicle moves can cross a coordinate-wise
 // local minimum without accepting an intermediate regional regression.
 if(!contactOK(best.m)){
  for(const step of [3,1.5,.75]){
   const from=a.slice();
   for(const [i,j] of [[0,4],[1,3],[2,3],[0,1],[1,2]])for(const si of [-1,1])for(const sj of [-1,1]){
    const c=from.slice();c[i]=clamp(c[i]+si*step,-limits[i],limits[i]);c[j]=clamp(c[j]+sj*step,-limits[j],limits[j]);
    const t=evaluate(c);if(t.total<best.total-1e-7&&t.raw<best.raw-1e-7){a=c;best=t;}
   }
  }
 }
 const improved=Number.isFinite(best.raw)&&best.raw<beforeScore-1e-7&&(
  (best.m.missing||0)<(beforeM.missing||0) || best.raw<beforeScore*.985 ||
  maxOut(best.m)<maxOut(beforeM)-.75 ||
  lowerOut(best.m)<lowerOut(beforeM)-.75
 );
 if(!improved){
  return {...input,v72:{active:false,reason:'no-safe-improvement',before:beforeM,after:beforeM,beforeScore,afterScore:beforeScore,params:[0,0,0,0,0,0],correctionDeg:0}};
 }
 const correctionDeg=v67?.qangle?v67.qangle(input.qS,best.p.qS):Math.hypot(a[0],a[1],a[2]);
 return {...best.p,
  before:input.before,after:input.after,
  v66Before:input.v66Before,v66After:input.v66After,v67:input.v67||null,
  fallback:input.fallback||null,
  v72:{active:true,reason:contactOK(best.m)?'contact-restored':'contact-improved',
   before:beforeM,after:best.m,beforeScore,afterScore:best.raw,params:a,correctionDeg,audit:best.audit||null}
 };
}

root.ScapulothoracicContactGuardV072=Object.freeze({
 version:'0.72.0',CFG,violation,noRegionalRegression,contactScore,contactOK,hardKinematicOK,refine
});
})(globalThis);
