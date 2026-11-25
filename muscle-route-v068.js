/* MuscleRouteEngineV068
 * Obstacle-aware refinement of v0.64 regional muscle guide curves.
 * Uses rest-relative thorax clearance and humeral-head keepout proxies.
 * It is not a collision/contact, force, fascia, tendon or finite-element solver.
 */
(function(root){
'use strict';
const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
const smooth=x=>{x=clamp(x,0,1);return x*x*(3-2*x)};
const add=(a,b)=>a.map((x,i)=>x+b[i]),sub=(a,b)=>a.map((x,i)=>x-b[i]),mul=(a,s)=>a.map(x=>x*s);
const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0),len=a=>Math.hypot(...a),norm=a=>mul(a,1/Math.max(len(a),1e-12));
const mix=(a,b,t)=>add(mul(a,1-t),mul(b,t));

const CFG=Object.freeze({
 v58RhomboidMajor:{thorax:[0,.92],inward:6,outward:14,shift:20},
 v58RhomboidMinor:{thorax:[0,.92],inward:6,outward:14,shift:20},
 v58Levator:{thorax:[0,.72],inward:7,outward:16,shift:22},
 v58Subscap:{head:[.28,1],headOut:12,shift:18},
 v59TrapUpper:{thorax:[0,.78],inward:7,outward:18,shift:24},
 v59TrapMiddle:{thorax:[0,.94],inward:6,outward:15,shift:22},
 v59TrapLower:{thorax:[0,.94],inward:6,outward:15,shift:22},
 v59Lat:{thorax:[0,.66],head:[.45,1],inward:7,outward:18,headOut:15,shift:26},
 v60PecClav:{head:[.42,1],headOut:14,shift:22},
 v60PecStern:{thorax:[0,.58],head:[.45,1],inward:6,outward:16,headOut:14,shift:24},
 v60PecAbd:{thorax:[0,.64],head:[.42,1],inward:7,outward:18,headOut:15,shift:25},
 v60PecMinor:{thorax:[0,.78],inward:5,outward:12,shift:18},
 v60Subclavius:{thorax:[0,.55],inward:5,outward:10,shift:14},
 v60DeltAnt:{head:[.08,.94],headOut:13,shift:21},
 v60DeltMid:{head:[.06,.96],headOut:13,shift:21},
 v60DeltPost:{head:[.08,.94],headOut:13,shift:21}
});

function mm(v,scale){return v*scale}
function active(range,t){return range&&t>=range[0]&&t<=range[1]}
function capShift(rest,p,original,max){
 const d=sub(p,original),n=len(d);return n<=max?p:add(original,mul(d,max/Math.max(n,1e-12)));
}
function thoraxCorrect(rest,p,cfg,env,t){
 if(!active(cfg.thorax,t)||!env.gapFn||!env.normalFn)return p;
 const g0=env.gapFn(rest),g1=env.gapFn(p),n=env.normalFn(p);
 if(!Number.isFinite(g0)||!Number.isFinite(g1)||!n||!n.every(Number.isFinite))return p;
 const lo=g0-mm(cfg.inward??7,env.scale),hi=g0+mm(cfg.outward??16,env.scale),target=clamp(g1,lo,hi);
 return add(p,mul(norm(n),target-g1));
}
function headCorrect(rest,p,cfg,env,t){
 if(!active(cfg.head,t)||!env.headRest||!env.headNow||!Number.isFinite(env.headRadius))return p;
 const r0=len(sub(rest,env.headRest)),v=sub(p,env.headNow),r1=len(v);if(r1<1e-9)return p;
 const clearance=mm(1.5,env.scale),minR=Math.min(r0,env.headRadius+clearance),maxR=r0+mm(cfg.headOut??14,env.scale);
 const target=clamp(r1,minR,maxR);return add(env.headNow,mul(v,target/r1));
}
function correction(rest,p,cfg,env,t,original=p){
 let q=p;
 // Order matters near the axilla: broad thorax support first, local humeral
 // head keepout second, as in ordered wrap-object geometry-path systems.
 q=thoraxCorrect(rest,q,cfg,env,t);
 q=headCorrect(rest,q,cfg,env,t);
 return capShift(rest,q,original,mm(cfg.shift??22,env.scale));
}
function violation(rest,p,cfg,env,t){
 let v=0,n=0;
 if(active(cfg.thorax,t)&&env.gapFn){
  const a=env.gapFn(rest),b=env.gapFn(p);
  if(Number.isFinite(a)&&Number.isFinite(b)){
   const lo=a-mm(cfg.inward??7,env.scale),hi=a+mm(cfg.outward??16,env.scale);
   const d=b<lo?lo-b:b>hi?b-hi:0;v+=(d/env.scale)**2;n++;
  }
 }
 if(active(cfg.head,t)&&env.headRest&&env.headNow&&Number.isFinite(env.headRadius)){
  const a=len(sub(rest,env.headRest)),b=len(sub(p,env.headNow));
  const lo=Math.min(a,env.headRadius+mm(1.5,env.scale)),hi=a+mm(cfg.headOut??14,env.scale);
  const d=b<lo?lo-b:b>hi?b-hi:0;v+=(d/env.scale)**2;n++;
 }
 return n?v/n:0;
}
function refine(id,field,posed,env){
 const cfg=CFG[id];if(!cfg||!field?.controls||!posed||!env?.scale)return {posed,stats:{enabled:false,id,reason:'unregistered'}};
 const base=posed.map(row=>row.map(p=>p.slice())),out=posed.map(row=>row.map(p=>p.slice()));
 let before=0,after=0,count=0,maxShift=0;
 // Correct only the two interior Bezier stations. Attachments remain exact.
 for(let lane=0;lane<out.length;lane++)for(let j=1;j<=2;j++){
  const t=j/3,rest=field.controls[lane][j],p=base[lane][j];
  before+=violation(rest,p,cfg,env,t);count++;
  out[lane][j]=correction(rest,p,cfg,env,t,p);
 }
 // Neighbor smoothing suppresses lane-to-lane ridges without changing end attachments.
 const tmp=out.map(row=>row.map(p=>p.slice()));
 for(let lane=0;lane<out.length;lane++)for(let j=1;j<=2;j++){
  const neighbors=[tmp[lane][j]];if(lane>0)neighbors.push(tmp[lane-1][j]);if(lane+1<tmp.length)neighbors.push(tmp[lane+1][j]);
  let avg=[0,0,0];for(const p of neighbors)avg=add(avg,p);avg=mul(avg,1/neighbors.length);
  let q=mix(tmp[lane][j],avg,.18);
  q=correction(field.controls[lane][j],q,cfg,env,j/3,base[lane][j]);
  out[lane][j]=q;
 }
 for(let lane=0;lane<out.length;lane++)for(let j=1;j<=2;j++){
  after+=violation(field.controls[lane][j],out[lane][j],cfg,env,j/3);count++;
  maxShift=Math.max(maxShift,len(sub(out[lane][j],base[lane][j]))/env.scale);
 }
 before/=Math.max(1,count/2);after/=Math.max(1,count/2);
 // Never worsen obstacle score; leave v0.64/v0.65 deformation untouched otherwise.
 if(!Number.isFinite(after)||after>before+1e-9)return {posed:base,stats:{enabled:true,id,adopted:false,before,after:before,maxShiftMm:0}};
 return {posed:out,stats:{enabled:true,id,adopted:after<before-1e-9,before,after,maxShiftMm:maxShift}};
}
function audit(id,field,posed,env){return refine(id,field,posed,env).stats}
root.MuscleRouteEngineV068=Object.freeze({version:'0.68.0',config:CFG,refine,audit});
})(globalThis);
