import '../serratus-v065.js';
const E=globalThis.SerratusEngineV065,assert=(x,m)=>{if(!x)throw Error(m)};
const verts=[],faces=[],weights=[],ribs=[];
for(let k=0;k<9;k++)ribs.push({level:k+1,points:Array.from({length:20},(_,j)=>[.02*j,1.7-k*.12,.38+.02*Math.sin(j)])});
for(let i=0;i<=18;i++)for(let j=0;j<18;j++){
 const t=i/18,s=j/17*8;verts.push([.2+.55*t,1.7-s*.12,.42-.35*t+.025*Math.sin(s)]);weights.push(t);
 if(i<18&&j<17){const k=i*18+j;faces.push(k,k+1,k+18,k+1,k+19,k+18)}
}
const medial=[[.78,1.69,.08],[.76,1.5,.06],[.72,1.05,.04],[.66,.72,.03]];
const F=E.build(verts,faces,weights,ribs,medial);assert(F&&F.controls.length===9,'9-slip build');
for(let i=1;i<F.slip.length;i++)assert(Number.isFinite(F.slip[i])&&F.slip[i]>=0&&F.slip[i]<=8,'finite slip coordinate');
const A={q:[0,0,0,1],p:[0,0,0],rest:[0,0,0]},P=E.pose(F,A,A);let neutral=0;
for(let i=0;i<verts.length;i++){const q=E.point(verts[i],weights[i],F.slip[i],F,P,A,A,verts[i]);neutral=Math.max(neutral,Math.hypot(...q.map((x,j)=>x-verts[i][j])))}
assert(neutral<1e-8,'neutral preservation');
let poses=0,adopted=0;
for(const deg of [30,60,90,120,150]){
 const a=deg*Math.PI/360,b={q:[0,0,Math.sin(a),Math.cos(a)],p:[.78,1.69,.08],rest:[.78,1.69,.08]};
 const base=(p,t)=>{const d=p.map((x,k)=>x-b.rest[k]),c=Math.cos(2*a),s=Math.sin(2*a),q=[b.p[0]+c*d[0]-s*d[1],b.p[1]+s*d[0]+c*d[1],b.p[2]+d[2]];return p.map((x,k)=>x*(1-t)+q[k]*t)};
 const posed=E.pose(F,A,b),r=E.assess(F,posed,A,b,base),r2=E.assess(F,posed,A,b,base);
 assert(r.score<=r.baseline+1e-9,'non-worsening gate');assert(JSON.stringify(r)===JSON.stringify(r2),'determinism');
 if(r.alpha>0)adopted++;poses++;
}
console.log({version:E.version,slips:F.controls.length,edges:F.refs.length,poses,adopted,neutralError:neutral});
