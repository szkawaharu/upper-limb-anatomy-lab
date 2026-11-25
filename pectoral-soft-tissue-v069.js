/* PectoralSoftTissueEngineV069
 * Pectoralis-major-specific soft-tissue restraint layered after v0.64/v0.68.2.
 * Not a skin, fascia, contact-pressure, tendon-force or FEM solver.
 */
(function(root){
'use strict';
const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
const smooth=x=>{x=clamp(x,0,1);return x*x*(3-2*x)};
const add=(a,b)=>a.map((x,i)=>x+b[i]);
const sub=(a,b)=>a.map((x,i)=>x-b[i]);
const mul=(a,s)=>a.map(x=>x*s);
const len=a=>Math.hypot(...a);
const norm=a=>mul(a,1/Math.max(len(a),1e-12));
const mix=(a,b,t)=>add(mul(a,1-t),mul(b,t));

const CFG=Object.freeze({
 v60PecClav:{fasciaFull:.46,fasciaEnd:.76,outMm:8,inMm:5,maxFasciaShiftMm:8,tendonStart:.70,baselineBlend:.62,bridgeStart:.48,bridgeEnd:.93,bridgeMm:9},
 v60PecStern:{fasciaFull:.50,fasciaEnd:.76,outMm:9,inMm:5,maxFasciaShiftMm:9,tendonStart:.66,baselineBlend:.66,bridgeStart:.46,bridgeEnd:.94,bridgeMm:10},
 v60PecAbd:{fasciaFull:.52,fasciaEnd:.78,outMm:10,inMm:6,maxFasciaShiftMm:10,tendonStart:.63,baselineBlend:.68,bridgeStart:.44,bridgeEnd:.95,bridgeMm:11}
});
function mm(x,scale){return x*scale}
function curve(c,t){const u=1-t;return add(add(mul(c[0],u*u*u),mul(c[1],3*u*u*t)),add(mul(c[2],3*u*t*t),mul(c[3],t*t*t)))}
function laneCurve(controls,l){const n=controls.length;if(n===1)return controls[0];const x=clamp(l,0,n-1),i=Math.min(n-2,Math.floor(x)),f=x-i;return controls[i].map((p,j)=>mix(p,controls[i+1][j],f))}
function fasciaWeight(cfg,t){if(t<=cfg.fasciaFull)return 1;if(t>=cfg.fasciaEnd)return 0;return 1-smooth((t-cfg.fasciaFull)/(cfg.fasciaEnd-cfg.fasciaFull))}
function tendonWeight(cfg,t){return t<=cfg.tendonStart?0:smooth((t-cfg.tendonStart)/(1-cfg.tendonStart))}
function bridgeWeight(cfg,t){if(t<=cfg.bridgeStart||t>=cfg.bridgeEnd)return 0;const a=smooth((t-cfg.bridgeStart)/Math.max(.08,(cfg.bridgeEnd-cfg.bridgeStart)*.35)),b=1-smooth((t-(cfg.bridgeEnd-.12))/.12);return clamp(a*b,0,1)}
function capVector(v,maxLen){const n=len(v);return n<=maxLen?v:mul(v,maxLen/Math.max(n,1e-12))}
function applyFascia(rest,q,cfg,env,t,diag){
 const fw=fasciaWeight(cfg,t);if(fw<=0||!env.gapFn||!env.normalFn)return q;
 const g0=env.gapFn(rest),g1=env.gapFn(q),n=env.normalFn(q);
 if(!Number.isFinite(g0)||!Number.isFinite(g1)||!n||!n.every(Number.isFinite))return q;
 const lo=g0-mm(cfg.inMm,env.scale),hi=g0+mm(cfg.outMm,env.scale),target=clamp(g1,lo,hi),excess=Math.abs(g1-target)/env.scale;
 if(excess<=.75)return q;
 const ramp=.74*smooth((excess-.75)/5.5)*fw,raw=mul(norm(n),(target-g1)*ramp),d=capVector(raw,mm(cfg.maxFasciaShiftMm,env.scale));
 diag.fasciaShiftMm=Math.max(diag.fasciaShiftMm,len(d)/env.scale);diag.fasciaActive=Math.max(diag.fasciaActive,ramp);
 return add(q,d);
}
function applyBridge(base,q,cfg,env,t,diag){
 const w=bridgeWeight(cfg,t);if(w<=0)return q;
 const d=sub(q,base),r=len(d),limit=mm(cfg.bridgeMm,env.scale);if(r<=limit)return q;
 const excess=(r-limit)/env.scale,k=.58*w*smooth(excess/8),target=add(base,mul(d,limit/Math.max(r,1e-12))),out=mix(q,target,k);
 diag.bridgeShiftMm=Math.max(diag.bridgeShiftMm,len(sub(out,q))/env.scale);return out;
}
function applyTendon(rest,q,base,cfg,env,t,l,field,posed,diag){
 const tw=tendonWeight(cfg,t);if(tw<=0)return q;
 const restC=curve(laneCurve(field.controls,l),t),nowC=curve(laneCurve(posed,l),t);
 let out=mix(q,base,cfg.baselineBlend*tw);
 const r0=len(sub(rest,restC)),v=sub(out,nowC),r=len(v);
 if(r>1e-10){const maxR=r0*(1+.08*(1-tw))+mm(.8,env.scale);if(r>maxR){const target=add(nowC,mul(v,maxR/r));out=mix(out,target,.82*tw)}}
 diag.tendonBlend=Math.max(diag.tendonBlend,tw);diag.tendonShiftMm=Math.max(diag.tendonShiftMm,len(sub(out,q))/env.scale);return out;
}
function refinePoint(id,rest,candidate,base,t,l,field,posed,env){
 const cfg=CFG[id];
 if(!cfg||!env?.scale||!field?.controls||!posed)return {point:candidate,diag:{enabled:false,id}};
 if(t<=0||t>=1)return {point:base,diag:{enabled:true,id,fasciaShiftMm:0,bridgeShiftMm:0,tendonShiftMm:0,fasciaActive:0,tendonBlend:t>=1?1:0,totalShiftMm:len(sub(base,candidate))/env.scale}};
 const diag={enabled:true,id,fasciaShiftMm:0,bridgeShiftMm:0,tendonShiftMm:0,fasciaActive:0,tendonBlend:0,totalShiftMm:0};
 let q=candidate;q=applyFascia(rest,q,cfg,env,t,diag);q=applyBridge(base,q,cfg,env,t,diag);q=applyTendon(rest,q,base,cfg,env,t,l,field,posed,diag);
 const maxTotal=mm(14+5*(1-t),env.scale),dv=sub(q,base);if(len(dv)>maxTotal)q=add(base,capVector(dv,maxTotal));
 diag.totalShiftMm=len(sub(q,candidate))/env.scale;return {point:q,diag};
}
function audit(id,field,posed,samples,env,pointFn){
 if(!CFG[id]||!field?.controls||!posed||!Array.isArray(samples)||typeof pointFn!=='function')return {enabled:false,id};
 let n=0,maxFascia=0,maxBridge=0,maxTendon=0,maxTotal=0,sumTotal=0,activeFascia=0,activeTendon=0;
 for(const s of samples){
  const base=pointFn(s,'base'),candidate=pointFn(s,'candidate'),r=refinePoint(id,s.p,candidate,base,s.w,s.l,field,posed,env),d=r.diag;if(!d?.enabled)continue;
  n++;maxFascia=Math.max(maxFascia,d.fasciaShiftMm||0);maxBridge=Math.max(maxBridge,d.bridgeShiftMm||0);maxTendon=Math.max(maxTendon,d.tendonShiftMm||0);maxTotal=Math.max(maxTotal,d.totalShiftMm||0);sumTotal+=d.totalShiftMm||0;
  if((d.fasciaActive||0)>.02)activeFascia++;if((d.tendonBlend||0)>.02)activeTendon++;
 }
 return {enabled:true,id,n,maxFasciaShiftMm:maxFascia,maxBridgeShiftMm:maxBridge,maxTendonShiftMm:maxTendon,maxTotalShiftMm:maxTotal,meanTotalShiftMm:n?sumTotal/n:0,fasciaActiveSamples:activeFascia,tendonActiveSamples:activeTendon};
}
root.PectoralSoftTissueEngineV069=Object.freeze({version:'0.69.0',config:CFG,refinePoint,audit,fasciaWeight,tendonWeight,bridgeWeight,curve,laneCurve});
})(globalThis);
