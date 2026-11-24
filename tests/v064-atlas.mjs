// Usage: node tests/v064-atlas.mjs /path/to/BodyParts3D/assets/BodyParts3D_data/stl
// Use source commit 2c83c698788cff57f6317bddcf712e7acfbea7d8.
// Simplified landmark weights and humeral-only rotation; not production validation.
import '../muscle-fiber-v064.js';
import {readFileSync} from 'node:fs';
import {join} from 'node:path';
const E=globalThis.FiberEngineV064,dir=process.argv[2];
if(!dir)throw Error('Supply directory containing the pinned source STL files');
const specs=[["FMA45874.stl",377384,"大胸筋腹部",[[-20,-155,1150],[-55,-155,1080]],[[-184,-97,1311],[-189,-97,1274]]],["FMA34690.stl",466584,"大胸筋鎖骨部",[[-18,-86,1343],[-80,-91,1344]],[[-184,-97,1311],[-189,-97,1274]]]];
const atlas=p=>[-(p[0]+14)*.00555,1.7+(p[2]-1341)*.00555,(p[1]+72)*.00555+.08],sub=(a,b)=>a.map((v,k)=>v-b[k]),dot=(a,b)=>a.reduce((s,v,k)=>s+v*b[k],0),dist=(p,a,b)=>{let v=sub(b,a),t=Math.max(0,Math.min(1,dot(sub(p,a),v)/dot(v,v)));return Math.hypot(...p.map((x,k)=>x-a[k]-v[k]*t))};
const results=[];
for(const [file,bytes,name,O,I] of specs){
 const buf=readFileSync(join(dir,file)),d=new DataView(buf.buffer,buf.byteOffset,buf.byteLength),count=d.getUint32(80,true);
 if(buf.length!==bytes||bytes!==84+50*count)throw Error('STL length '+file);
 const verts=[],faces=[],map=new Map(),w=[],o=O.map(atlas),i=I.map(atlas);
 for(let k=0;k<count;k++){let ids=[];for(let j=0;j<3;j++){let z=96+50*k+12*j,p=atlas([d.getFloat32(z,true),d.getFloat32(z+4,true),d.getFloat32(z+8,true)]),key=p.map(x=>Math.round(x/.00555)).join(",");if(!map.has(key)){map.set(key,verts.length);verts.push(p);let a=dist(p,...o),b=dist(p,...i),t=a<.00555*8&&a<b*.7?0:b<.00555*8&&b<a*.7?1:Math.max(0,Math.min(1,(a/Math.max(a+b,1e-9)-.10)/.80));w.push(t*t*(3-2*t));}ids.push(map.get(key))}if(new Set(ids).size===3)faces.push(...ids)}
 const f=E.build(verts,faces,w,1);if(!f)throw Error('guide unavailable '+file);
 const A={q:[0,0,0,1],p:[0,0,0],rest:[0,0,0]},gh=atlas([-172,-72,1320]);let improved=0;const rows=[];
 for(const angle of [0,45,90,135,180]){
  const a=angle*Math.PI/360,B={q:[0,0,Math.sin(a),Math.cos(a)],p:gh,rest:gh};
  const base=(p,t)=>{const x=p[0]-gh[0],y=p[1]-gh[1],c=Math.cos(a*2),s=Math.sin(a*2),q=[gh[0]+c*x-s*y,gh[1]+s*x+c*y,p[2]];return p.map((x,k)=>x*(1-t)+q[k]*t)};
  const posed=E.pose(f,A,B),q=E.assess(f,posed,A,B,base);
  if(q.score>q.baseline+1e-9||q.areaScore>q.areaBaseline+1e-8)throw Error('gate failed');
  if(q.alpha>0)improved++;
  rows.push({angle,alpha:+q.alpha.toFixed(3),base:+q.baseline.toFixed(5),score:+q.score.toFixed(5),areaBase:+q.areaBaseline.toFixed(5),areaScore:+q.areaScore.toFixed(5)});
 }
 results.push({name,file,vertices:verts.length,improved,results:rows});
}
console.log(JSON.stringify(results,null,2));
