import '../muscle-route-v068.js';

const E=globalThis.MuscleRouteEngineV068;
const assert=(x,m)=>{if(!x)throw Error(m)};
const controls=Array.from({length:6},(_,i)=>[[0,0,i*.01],[.30,.03,i*.01],[.68,.02,i*.01],[1,0,i*.01]]);
const field={controls},scale=.00555,gapFn=p=>p[2],normalFn=p=>[0,0,1];
let cases=0,adopted=0,maxShiftMm=0;

for(const id of Object.keys(E.config))for(const push of [-.30,-.18,-.08,.08,.22]){
 const posed=controls.map(r=>r.map(p=>p.slice()));
 for(const row of posed){row[1][2]+=push;row[2][2]+=push*.85}
 const env={scale,gapFn,normalFn,headRest:[1,0,0],headNow:[1,.02,0],headRadius:.13};
 const out=E.refine(id,field,posed,env);
 assert(out.posed.flat(2).every(Number.isFinite),id+' finite');
 if(out.stats.enabled)assert(out.stats.after<=out.stats.before+1e-9,id+' non-worsening obstacle score');
 for(let lane=0;lane<6;lane++){
  assert(JSON.stringify(out.posed[lane][0])===JSON.stringify(posed[lane][0]),id+' origin fixed');
  assert(JSON.stringify(out.posed[lane][3])===JSON.stringify(posed[lane][3]),id+' insertion fixed');
 }
 assert((out.stats.maxShiftMm||0)<=(E.config[id].shift??22)+1e-6,id+' bounded guide shift');
 if(out.stats.adopted)adopted++;
 maxShiftMm=Math.max(maxShiftMm,out.stats.maxShiftMm||0);cases++;
}
console.log({version:E.version,registered:Object.keys(E.config).length,cases,adopted,maxShiftMm,contactValidated:false});
