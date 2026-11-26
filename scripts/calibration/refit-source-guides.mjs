// Positive centroid control polygons, fitted to the FINAL harmonic material
// field. Do not independently project t then smooth it away from the fitted
// geometry. Atlas vertices and the fixed attachment masks are unchanged.
import{readFileSync,writeFileSync,mkdirSync}from'node:fs';
import '../../muscle-fiber-v064.js';import '../../pectoral-fan-v070.js';
const input=process.argv[2]||'test-results/v076-source-attached/bindings.json',out=process.argv[3]||'test-results/v076-source-refit';mkdirSync(out,{recursive:true});const all=JSON.parse(readFileSync(input));
for(const B of all){const vertices=Array.from({length:B.weights.length},(_,i)=>B.rest.slice(i*3,i*3+3)),w=B.weights,lo=Math.min(...w.slice(0,60000)),hi=Math.max(...w.slice(0,60000));let min=1,max=0;for(const x of w){min=Math.min(min,x);max=Math.max(max,x)}const t=w.map(x=>(x-min)/(max-min));
 let f=B.sheet?PectoralFanEngineV070.build(vertices,B.faces,t):null;
 if(!f&&B.spec.id!=='v58Serratus')for(const axis of [1,0,2]){f=FiberEngineV064.build(vertices,B.faces,t,axis);if(f)break;}
 if(B.spec.id==='v58Serratus')f={controls:B.field.controls,uv:B.lane,u:B.reference};
 if(!f)throw Error('No source guide '+B.spec.id);
 const lanes=Array.from(f.uv||f.slip),count=f.controls.length,cells=Array.from({length:count},()=>Array.from({length:4},()=>({s:[0,0,0],w:0})));for(let i=0;i<w.length;i++){const l=lanes[i],u=t[i]*3,p=vertices[i];for(let j=0;j<4;j++){const wt=Math.max(0,1-Math.abs(u-j)/1.15)**2;if(!wt)continue;for(let k=Math.max(0,Math.floor(l)-1);k<=Math.min(count-1,Math.ceil(l)+1);k++){const a=wt*Math.max(0,1-Math.abs(l-k)/1.4)**2,c=cells[k][j];c.w+=a;for(let h=0;h<3;h++)c.s[h]+=a*p[h];}}}
 const controls=cells.map((row,k)=>row.map((cell,j)=>{if(cell.w)return cell.s.map(x=>x/cell.w);let nearest=-1;for(let n=0;n<count;n++)if(cells[n][j].w&&(nearest<0||Math.abs(n-k)<Math.abs(nearest-k)))nearest=n;if(nearest<0)throw Error('Empty source material station '+B.spec.id);const c=cells[nearest][j];return c.s.map(x=>x/c.w)}));
 B.field={controls};B.lane=lanes;B.material=t;B.reference=f.u||B.reference;B.blendRange=[min,max];B.fieldMethod='source-harmonic-positive-control-polygons';B.guidePreparation={method:B.fieldMethod,independentClosestProjection:false,controlsAreInterpolationFit:false,anatomicalFasciclesMeasured:false};console.log(B.spec.id,min,max,count);
}
writeFileSync(out+'/bindings.json',JSON.stringify(all));
