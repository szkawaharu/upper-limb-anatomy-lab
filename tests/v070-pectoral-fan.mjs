import '../pectoral-fan-v070.js';

const E=globalThis.PectoralFanEngineV070;
const assert=(x,m)=>{if(!x)throw Error(m)};
const add=(a,b)=>a.map((x,i)=>x+b[i]);
const sub=(a,b)=>a.map((x,i)=>x-b[i]);
const mul=(a,s)=>a.map(x=>x*s);
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const len=a=>Math.hypot(...a);
const qrot=(q,v)=>{const u=q.slice(0,3),t=mul(cross(u,v),2);return add(v,add(mul(t,q[3]),cross(u,t)))};
const qa=(axis,a)=>{const n=len(axis)||1,s=Math.sin(a/2);return [axis[0]/n*s,axis[1]/n*s,axis[2]/n*s,Math.cos(a/2)]};

const vertices=[],weights=[],faces=[],nu=11,nt=18;
for(let j=0;j<nt;j++){
 const t=j/(nt-1);
 for(let i=0;i<nu;i++){
  const u=i/(nu-1)*2-1;
  // Broad origin, narrow distal tendon, slight sheet camber/thickness.
  vertices.push([t,(1-t)*u*.58,.022*Math.sin(Math.PI*t)*(1-u*u)]);
  weights.push(t);
 }
}
for(let j=0;j<nt-1;j++)for(let i=0;i<nu-1;i++){
 const a=j*nu+i,b=a+1,c=a+nu,d=c+1;faces.push(a,c,b,b,c,d);
}

const field=E.build(vertices,faces,weights);
assert(field&&field.engine==='fan070','field builds');
assert(field.controls.length===6&&field.controls.every(r=>r.length===4),'6x4 sheet controls');
assert(field.uv.length===vertices.length,'per-vertex fan lane');

const originWidth=len(sub(E.curve(field.controls[5],.12),E.curve(field.controls[0],.12)));
const distalWidth=len(sub(E.curve(field.controls[5],.88),E.curve(field.controls[0],.88)));
assert(originWidth>distalWidth*1.8,'broad origin converges distally');

const ID=[0,0,0,1],A={q:ID,p:[0,0,0],rest:[0,0,0]};
let previous=null,maxAdjacentJump=0,minLaneSpacing=Infinity,maxNormalFlip=1,cases=0;
for(let k=0;k<=120;k++){
 const angle=1.65*k/120,B={q:qa([0,0,1],angle),p:[.14,.09,.03],rest:[1,0,0]};
 const posed=E.pose(field,A,B,null);
 const baseline=(p,t)=>add(mul(p,1-t),mul(add(B.p,qrot(B.q,sub(p,B.rest))),t));

 // Sheet lanes stay separated through the belly.
 for(const t of [.12,.25,.45,.65,.82])for(let lane=0;lane<5;lane++){
  const p=E.curve(posed[lane],t),q=E.curve(posed[lane+1],t);
  minLaneSpacing=Math.min(minLaneSpacing,len(sub(p,q)));
 }
 // Local frame normal must not randomly invert along a lane.
 for(let lane=0;lane<6;lane++){
  let prevN=null;
  for(const t of [.18,.3,.45,.6,.78]){
   const F=E.frame(posed,lane,t);assert(F,'valid sheet frame');
   if(prevN)maxNormalFlip=Math.min(maxNormalFlip,prevN[0]*F.n[0]+prevN[1]*F.n[1]+prevN[2]*F.n[2]);
   prevN=F.n;
  }
 }

 const now={};
 for(let i=0;i<vertices.length;i+=9){
  const p=E.point(vertices[i],weights[i],field.uv[i],field,posed,A,B,baseline(vertices[i],weights[i]),1);
  assert(p.every(Number.isFinite),'finite point');
  if(weights[i]===0||weights[i]===1)assert(JSON.stringify(p)===JSON.stringify(baseline(vertices[i],weights[i])),'attachments exact');
  now[i]=p;
  if(previous?.[i])maxAdjacentJump=Math.max(maxAdjacentJump,len(sub(p,previous[i])));
 }
 previous=now;cases++;
}
assert(minLaneSpacing>.01,'no fan lane collapse');
assert(maxAdjacentJump<.05,'continuous elevation sweep');
assert(maxNormalFlip>.2,'sheet normal continuity');

console.log({version:E.version,cases,originWidth,distalWidth,minLaneSpacing,maxAdjacentJump,maxNormalFlip});
