import '../pectoral-soft-tissue-v069.js';

const E=globalThis.PectoralSoftTissueEngineV069;
const assert=(x,m)=>{if(!x)throw Error(m)};
const add=(a,b)=>a.map((x,i)=>x+b[i]);
const sub=(a,b)=>a.map((x,i)=>x-b[i]);
const len=a=>Math.hypot(...a);

const field={controls:Array.from({length:6},(_,i)=>[
 [0,0,i*.01],[.28,.01,i*.01],[.68,.03,i*.01],[1,0,i*.01]
])};
const posed=field.controls.map(r=>r.map(p=>p.slice()));
for(const r of posed){r[1][2]+=.03;r[2][2]+=.08}
const env={scale:.00555,gapFn:p=>p[2],normalFn:p=>[0,0,1]};

let cases=0,maxAdjacentJumpMm=0,fasciaReductions=0,tendonReductions=0;
for(const id of Object.keys(E.config)){
 let previous=null;
 for(let k=0;k<=200;k++){
  const t=k/200,l=2.4;
  const rest=E.curve(E.laneCurve(field.controls,l),t);
  const currentCenter=E.curve(E.laneCurve(posed,l),t);
  const base=add(currentCenter,[0,.012,.010]);
  const balloon=add(currentCenter,[0,.020,.080*(Math.sin(Math.PI*t)**2)+.020]);
  const out=E.refinePoint(id,rest,balloon,base,t,l,field,posed,env);

  assert(out.point.every(Number.isFinite),id+' finite');
  if(t===0||t===1)assert(JSON.stringify(out.point)===JSON.stringify(base),id+' exact attachment');
  if(previous)maxAdjacentJumpMm=Math.max(maxAdjacentJumpMm,len(sub(out.point,previous))/env.scale);

  if(t<.55&&env.gapFn(out.point)<env.gapFn(balloon)-1e-10)fasciaReductions++;
  if(t>.82&&len(sub(out.point,base))<len(sub(balloon,base))-1e-10)tendonReductions++;
  previous=out.point;cases++;
 }
}
assert(maxAdjacentJumpMm<3,'continuous pose/tendon transition');
assert(fasciaReductions>0,'fascia restraint active');
assert(tendonReductions>0,'tendon preservation active');

const untouched=E.refinePoint('not-a-pectoralis',[0,0,0],[1,2,3],[0,0,0],.5,2,field,posed,env);
assert(JSON.stringify(untouched.point)===JSON.stringify([1,2,3]),'non-pectoralis untouched');

console.log({version:E.version,cases,maxAdjacentJumpMm,fasciaReductions,tendonReductions,registered:Object.keys(E.config).length});
