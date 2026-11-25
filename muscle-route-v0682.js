/* MuscleRouteEngineV0682
 * Conservative, continuous obstacle-aware refinement of v0.64 muscle guides.
 * v0.68.2 fixes two problems from v0.68:
 *  1) lane smoothing is applied to correction VECTORS, never to the original
 *     control-point positions, so a muscle with no obstacle violation is unchanged;
 *  2) obstacle correction fades in continuously with a deadband and is limited by
 *     control-polygon bend, reducing pose-to-pose shape snapping and S-bends.
 *
 * This remains an educational geometry-path approximation, not a contact/force/FEM solver.
 */
(function(root){
'use strict';
const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
const smooth=x=>{x=clamp(x,0,1);return x*x*(3-2*x)};
const add=(a,b)=>a.map((x,i)=>x+b[i]),sub=(a,b)=>a.map((x,i)=>x-b[i]),mul=(a,s)=>a.map(x=>x*s);
const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0),len=a=>Math.hypot(...a),norm=a=>mul(a,1/Math.max(len(a),1e-12));
const mix=(a,b,t)=>add(mul(a,1-t),mul(b,t));
const deg=r=>r*180/Math.PI;

const CFG=Object.freeze({
 v58RhomboidMajor:{thorax:[0,.92],inward:6,outward:14,shift:16},
 v58RhomboidMinor:{thorax:[0,.92],inward:6,outward:14,shift:16},
 v58Levator:{thorax:[0,.72],inward:7,outward:16,shift:17},
 v58Subscap:{head:[.28,1],headOut:12,shift:14},
 v59TrapUpper:{thorax:[0,.78],inward:7,outward:18,shift:17},
 v59TrapMiddle:{thorax:[0,.94],inward:6,outward:15,shift:16},
 v59TrapLower:{thorax:[0,.94],inward:6,outward:15,shift:16},
 v59Lat:{thorax:[0,.66],head:[.45,1],inward:7,outward:18,headOut:15,shift:18},
 v60PecClav:{head:[.42,1],headOut:14,shift:15},
 v60PecStern:{thorax:[0,.58],head:[.45,1],inward:6,outward:16,headOut:14,shift:17},
 v60PecAbd:{thorax:[0,.64],head:[.42,1],inward:7,outward:18,headOut:15,shift:17},
 v60PecMinor:{thorax:[0,.78],inward:5,outward:12,shift:14},
 v60Subclavius:{thorax:[0,.55],inward:5,outward:10,shift:11},
 v60DeltAnt:{head:[.08,.94],headOut:13,shift:15},
 v60DeltMid:{head:[.06,.96],headOut:13,shift:15},
 v60DeltPost:{head:[.08,.94],headOut:13,shift:15}
});

function mm(v,scale){return v*scale}
function active(range,t){return range&&t>=range[0]&&t<=range[1]}
function capDelta(d,max){const n=len(d);return n<=max?d:mul(d,max/Math.max(n,1e-12))}
function thoraxTarget(rest,p,cfg,env,t){
 if(!active(cfg.thorax,t)||!env.gapFn||!env.normalFn)return p;
 const g0=env.gapFn(rest),g1=env.gapFn(p),n=env.normalFn(p);
 if(!Number.isFinite(g0)||!Number.isFinite(g1)||!n||!n.every(Number.isFinite))return p;
 const lo=g0-mm(cfg.inward??7,env.scale),hi=g0+mm(cfg.outward??16,env.scale);
 return add(p,mul(norm(n),clamp(g1,lo,hi)-g1));
}
function headTarget(rest,p,cfg,env,t){
 if(!active(cfg.head,t)||!env.headRest||!env.headNow||!Number.isFinite(env.headRadius))return p;
 const r0=len(sub(rest,env.headRest)),v=sub(p,env.headNow),r1=len(v);if(r1<1e-9)return p;
 const clearance=mm(1.5,env.scale),minR=Math.min(r0,env.headRadius+clearance),maxR=r0+mm(cfg.headOut??14,env.scale);
 return add(env.headNow,mul(v,clamp(r1,minR,maxR)/r1));
}
function rawDelta(rest,p,cfg,env,t){
 let q=thoraxTarget(rest,p,cfg,env,t);
 q=headTarget(rest,q,cfg,env,t);
 return capDelta(sub(q,p),mm(cfg.shift??15,env.scale));
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
function bendDeg(row){
 const seg=[sub(row[1],row[0]),sub(row[2],row[1]),sub(row[3],row[2])];
 let m=0;
 for(let i=0;i<2;i++){
  const a=norm(seg[i]),b=norm(seg[i+1]);
  m=Math.max(m,deg(Math.acos(clamp(dot(a,b),-1,1))));
 }
 return m;
}
function curvatureLimit(baseRow,deltas,env){
 const baseB=bendDeg(baseRow),trial=baseRow.map((p,i)=>i===0||i===3?p:add(p,deltas[i]));
 const trialB=bendDeg(trial),limit=Math.max(34,baseB+16);
 if(!Number.isFinite(trialB)||trialB<=limit+1e-9)return 1;
 const excess=trialB-baseB;
 return clamp((limit-baseB)/Math.max(excess,1e-9),0,1);
}
function refine(id,field,posed,env){
 const cfg=CFG[id];
 if(!cfg||!field?.controls||!posed||!env?.scale)return {posed,stats:{enabled:false,id,reason:'unregistered'}};
 const base=posed.map(row=>row.map(p=>p.slice()));
 const deltas=base.map(row=>row.map(()=>[0,0,0]));
 let before=0,count=0,maxRawMm=0;

 for(let lane=0;lane<base.length;lane++)for(let j=1;j<=2;j++){
  const t=j/3,rest=field.controls[lane][j],p=base[lane][j],v=violation(rest,p,cfg,env,t);
  before+=v;count++;
  const severity=Math.sqrt(Math.max(0,v)); // mm because violation is normalized by scale
  // 1 mm deadband, then a soft 6 mm ramp. Maximum correction strength 68%.
  const activation=.68*smooth((severity-1)/6);
  const d=mul(rawDelta(rest,p,cfg,env,t),activation);
  deltas[lane][j]=d;maxRawMm=Math.max(maxRawMm,len(d)/env.scale);
 }
 before/=Math.max(1,count);

 // IMPORTANT: smooth only the displacement field. Never average the original
 // control points: a zero-correction muscle must remain byte-for-byte unchanged.
 const raw=deltas.map(row=>row.map(d=>d.slice()));
 for(let lane=0;lane<deltas.length;lane++)for(let j=1;j<=2;j++){
  const neighbors=[raw[lane][j]];if(lane>0)neighbors.push(raw[lane-1][j]);if(lane+1<raw.length)neighbors.push(raw[lane+1][j]);
  let avg=[0,0,0];for(const d of neighbors)avg=add(avg,d);avg=mul(avg,1/neighbors.length);
  deltas[lane][j]=mix(raw[lane][j],avg,.16);
 }

 // Limit new curvature lane-by-lane. Scaling the pair of interior correction
 // vectors together preserves the existing endpoint tangents as much as possible.
 let curveScaleMin=1;
 for(let lane=0;lane<base.length;lane++){
  const k=curvatureLimit(base[lane],deltas[lane],env);curveScaleMin=Math.min(curveScaleMin,k);
  deltas[lane][1]=mul(deltas[lane][1],k);deltas[lane][2]=mul(deltas[lane][2],k);
 }

 const candidate=base.map((row,lane)=>row.map((p,j)=>j===0||j===3?p.slice():add(p,deltas[lane][j])));
 let after=0,maxShift=0,maxBendIncrease=0;
 for(let lane=0;lane<candidate.length;lane++){
  maxBendIncrease=Math.max(maxBendIncrease,bendDeg(candidate[lane])-bendDeg(base[lane]));
  for(let j=1;j<=2;j++){
   after+=violation(field.controls[lane][j],candidate[lane][j],cfg,env,j/3);
   maxShift=Math.max(maxShift,len(sub(candidate[lane][j],base[lane][j]))/env.scale);
  }
 }
 after/=Math.max(1,count);

 // No meaningful violation => no geometry changes at all.
 if(before<=1+1e-9||maxShift<.05)
  return {posed:base,stats:{enabled:true,id,adopted:false,before,after:before,maxShiftMm:0,blend:0,maxBendIncreaseDeg:0}};

 // Continuous non-worsening blend. In the unusual case the coupled correction
 // worsens the score, reduce the full displacement uniformly instead of toggling
 // immediately between two unrelated guide shapes.
 let blend=1,output=candidate,outScore=after;
 if(!Number.isFinite(after)||after>before){
  let best=base,bestScore=before,bestBlend=0;
  for(const k of [.8,.6,.4,.2]){
   const trial=base.map((row,lane)=>row.map((p,j)=>j===0||j===3?p.slice():add(p,mul(deltas[lane][j],k))));
   let score=0;
   for(let lane=0;lane<trial.length;lane++)for(let j=1;j<=2;j++)score+=violation(field.controls[lane][j],trial[lane][j],cfg,env,j/3);
   score/=Math.max(1,count);
   if(score<bestScore){best=trial;bestScore=score;bestBlend=k}
  }
  output=best;outScore=bestScore;blend=bestBlend;
 }
 return {posed:output,stats:{enabled:true,id,adopted:blend>0&&outScore<before-1e-6,before,after:outScore,maxShiftMm:maxShift*blend,blend,maxBendIncreaseDeg:maxBendIncrease*blend,curveScaleMin,maxRawMm}};
}
function audit(id,field,posed,env){return refine(id,field,posed,env).stats}
root.MuscleRouteEngineV0682=Object.freeze({version:'0.68.2',config:CFG,refine,audit,bendDeg});
})(globalThis);
