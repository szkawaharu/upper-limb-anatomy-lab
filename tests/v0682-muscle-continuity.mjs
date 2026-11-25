import '../muscle-route-v0682.js';

const E=globalThis.MuscleRouteEngineV0682;
const assert=(x,m)=>{if(!x)throw Error(m)};
const controls=Array.from({length:6},(_,i)=>[[0,0,i*.01],[.30,.03,i*.01],[.68,.02,i*.01],[1,0,i*.01]]);
const field={controls},scale=.00555,gapFn=p=>p[2],normalFn=p=>[0,0,1];
const dist=(a,b)=>Math.hypot(a[0]-b[0],a[1]-b[1],a[2]-b[2]);
let cases=0,maxAdjacentJumpMm=0,maxBendIncreaseDeg=0,noViolationExact=0;

for(const id of Object.keys(E.config)){
 let previous=null;
 for(let n=0;n<=100;n++){
  const push=-.30+n*.0052;
  const posed=controls.map(r=>r.map(p=>p.slice()));
  for(const row of posed){row[1][2]+=push;row[2][2]+=push*.85}
  const env={scale,gapFn,normalFn,headRest:[1,0,0],headNow:[1,.02,0],headRadius:.13};
  const out=E.refine(id,field,posed,env);
  assert(out.posed.flat(2).every(Number.isFinite),id+' finite');
  assert(out.stats.after<=out.stats.before+1e-8,id+' non-worsening obstacle score');
  for(let lane=0;lane<6;lane++){
   assert(JSON.stringify(out.posed[lane][0])===JSON.stringify(posed[lane][0]),id+' origin fixed');
   assert(JSON.stringify(out.posed[lane][3])===JSON.stringify(posed[lane][3]),id+' insertion fixed');
  }
  if(out.stats.before<=1.000001){
   assert(JSON.stringify(out.posed)===JSON.stringify(posed),id+' zero/small violation must preserve original guide exactly');
   noViolationExact++;
  }
  maxBendIncreaseDeg=Math.max(maxBendIncreaseDeg,out.stats.maxBendIncreaseDeg||0);
  if(previous){
   for(let lane=0;lane<6;lane++)for(let j=1;j<=2;j++)
    maxAdjacentJumpMm=Math.max(maxAdjacentJumpMm,dist(out.posed[lane][j],previous[lane][j])/scale);
  }
  previous=out.posed;cases++;
 }
}
assert(maxAdjacentJumpMm<2,'continuous pose sweep: no guide snapping');
assert(maxBendIncreaseDeg<=16.01,'curvature limiter');
console.log({version:E.version,cases,noViolationExact,maxAdjacentJumpMm,maxBendIncreaseDeg});
