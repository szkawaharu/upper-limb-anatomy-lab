// Candidate ONLY: C1 interpolation of screened 15-degree node coefficients.
 // Smooth coefficient does not certify collision-free intermediate poses.
import {readFileSync} from 'node:fs';
export function createBlend(report){
 const planes=[0,40,90],nodes=Array.from({length:13},(_,i)=>15*i);
 const table=planes.map(plane=>{
  const part=report.results.find(r=>r.rows[0].plane===plane);
  if(!part||part.summary.poses!==13||part.summary.controlRejected.length)throw Error('Invalid coarse evidence');
  return nodes.map(elev=>{
   const a=part.summary.improved.find(r=>r.elev===elev);
   // Negligible numerical differences are not treated as positive evidence.
   return a&&a.controlRms-a.candidateRms>1e-6?a.alpha:0;
  });
 });
 const bracket=(values,x)=>{
  if(!Number.isFinite(x)||x<values[0]||x>values.at(-1))throw Error('Outside screened domain');
  const i=Math.min(values.length-2,values.findIndex((v,j)=>j>0&&x<=v)-1);
  if(i<0)throw Error('Invalid interpolation bracket');
  const t=(x-values[i])/(values[i+1]-values[i]);
  return [i,t*t*(3-2*t)];
 };
 const sample=(elev,plane)=>{
  const [i,u]=bracket(nodes,elev),[j,v]=bracket(planes,plane);
  const a=table[j][i]*(1-u)+table[j][i+1]*u,b=table[j+1][i]*(1-u)+table[j+1][i+1]*u;
  return a*(1-v)+b*v;
 };
 return {sample,table,nodes,planes,sourceRun:report.runId,productAccepted:false};
}
export function readBlend(path='docs/validation/v076/shoulder-range-screening.json'){
 return createBlend(JSON.parse(readFileSync(path)));
}
