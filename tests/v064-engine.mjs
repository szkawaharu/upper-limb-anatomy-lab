import '../muscle-fiber-v064.js';
const E=globalThis.FiberEngineV064,verts=[],faces=[],weights=[];
const assert=(v,msg)=>{if(!v)throw Error(msg)};
for(let i=0;i<=20;i++)for(let j=0;j<16;j++){
 const t=i/20,a=j*Math.PI/8;verts.push([t,.12*Math.cos(a),.12*Math.sin(a)]);weights.push(t);
 if(i<20){const k=i*16+j,n=i*16+(j+1)%16;faces.push(k,n,k+16,n,n+16,k+16)}
}
const F=E.build(verts,faces,weights,1);assert(F&&F.controls.length===6,'field');
const A={q:[0,0,0,1],p:[0,0,0],rest:[0,0,0]},P=E.pose(F,A,A);let neutral=0;
for(let i=0;i<verts.length;i++){const q=E.point(verts[i],weights[i],F.uv[i],F,P,A,A,verts[i]);neutral=Math.max(neutral,Math.hypot(...q.map((x,j)=>x-verts[i][j])))}
assert(neutral<1e-8,'neutral');
let count=0,adopted=0;
for(const d of [0,30,60,90,120,150,180]){
 const a=d*Math.PI/360,b={q:[0,0,Math.sin(a),Math.cos(a)],rest:[1,0,0],p:[1+.25*Math.sin(a),.3*Math.sin(a),0]};
 const base=(p,t)=>{const x=p[0]-1,y=p[1],c=Math.cos(2*a),s=Math.sin(2*a),q=[b.p[0]+c*x-s*y,b.p[1]+s*x+c*y,p[2]];return p.map((v,j)=>v*(1-t)+q[j]*t)};
 const pose=E.pose(F,A,b),r=E.assess(F,pose,A,b,base),r2=E.assess(F,pose,A,b,base);
 assert(r.score<=r.baseline+1e-9&&r.areaScore<=r.areaBaseline+1e-8,'edge/area gate');
 assert(JSON.stringify(r)===JSON.stringify(r2),'determinism');
 const constrained=E.assess(F,pose,A,b,base,(rest,p)=>Math.max(0,Math.abs(p[2])-.12));
 assert(constrained.contactScore<=constrained.contactBaseline+1e-8,'clearance gate');
 for(let i=0;i<verts.length;i++){
  const q=E.point(verts[i],weights[i],F.uv[i],F,pose,A,b,base(verts[i],weights[i]),r.alpha);
  assert(q.every(Number.isFinite),'finite');
  if(weights[i]===0||weights[i]===1)assert(Math.hypot(...q.map((v,j)=>v-base(verts[i],weights[i])[j]))<1e-10,'attachment');
 }
 if(r.alpha>0)adopted++;count++;
}
const a=[[0,0,0,.1],[.5,0,0,.2],[1,0,0,.1]],b=[[0,0,0,.1],[.6,0,0,.2],[1.2,0,0,.1]],corrected=E.tubeCorrect(b,E.tubeVolume(a));
assert(Math.abs(E.tubeVolume(corrected)/E.tubeVolume(a)-1)<.002,'tube volume');
assert(corrected[0][3]===b[0][3]&&corrected[2][3]===b[2][3],'tube endpoints');
assert(E.build(verts,faces,weights,7)===null,'unknown axis');
console.log({syntheticPoses:count,adopted,neutralError:neutral,tubeVolumeRatio:E.tubeVolume(corrected)/E.tubeVolume(a),areaSamples:F.tris.length});
