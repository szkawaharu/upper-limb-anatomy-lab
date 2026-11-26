/* ScapulothoracicEngineV066
 * Rigid scapula constrained against a thoracic surface proxy in three regions.
 * Geometry-first educational/research model. No cartilage, bursa, contact force,
 * ligament, muscle-force, pathology or clinical-normality solver.
 */
(function(root){
'use strict';
const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
const smooth=x=>{x=clamp(x,0,1);return x*x*(3-2*x)};
const add=(a,b)=>a.map((x,i)=>x+b[i]),sub=(a,b)=>a.map((x,i)=>x-b[i]),mul=(a,s)=>a.map(x=>x*s);
const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0),cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const len=a=>Math.hypot(...a),norm=a=>mul(a,1/Math.max(len(a),1e-12));
const qa=(axis,a)=>{const n=norm(axis),s=Math.sin(a/2);return [n[0]*s,n[1]*s,n[2]*s,Math.cos(a/2)]};
const qm=(a,b)=>[a[3]*b[0]+a[0]*b[3]+a[1]*b[2]-a[2]*b[1],a[3]*b[1]-a[0]*b[2]+a[1]*b[3]+a[2]*b[0],a[3]*b[2]+a[0]*b[1]-a[1]*b[0]+a[2]*b[3],a[3]*b[3]-a[0]*b[0]-a[1]*b[1]-a[2]*b[2]];
const qrot=(q,v)=>{const u=q.slice(0,3),t=mul(cross(u,v),2);return add(v,add(mul(t,q[3]),cross(u,t)))};
const transform=(q,o,p)=>add(o,qrot(q,p));
const rad=d=>d*Math.PI/180;

function regionOf(p,lo,hi){
 const y=clamp((p[1]-lo[1])/Math.max(hi[1]-lo[1],1e-9),0,1);
 return y>.66?'upper':y>.34?'middle':'lower';
}
function build(points,restAC,profileNormal){
 if(!points?.length||!restAC||typeof profileNormal!=='function')return null;
 const lo=[Infinity,Infinity,Infinity],hi=[-Infinity,-Infinity,-Infinity];
 for(const v of points)for(let i=0;i<3;i++){lo[i]=Math.min(lo[i],v.p[i]);hi[i]=Math.max(hi[i],v.p[i])}
 const groups={upper:[],middle:[],lower:[]};
 for(const v of points){const region=v.region||regionOf(v.p,lo,hi);groups[region].push({...v,region})}
 if(Object.values(groups).some(x=>x.length<2))return null;
 const centroid={};
 for(const [name,arr] of Object.entries(groups)){
  let c=[0,0,0],w=0;for(const v of arr){const k=v.weight||1;c=add(c,mul(add(restAC,v.p),k));w+=k}centroid[name]=mul(c,1/Math.max(w,1e-9));
 }
 // A stable whole-blade costal normal: medial-border superior->inferior crossed
 // with medial->lateral direction. The sign is irrelevant; tangency uses abs(dot).
 const upper=centroid.upper,lower=centroid.lower,mid=centroid.middle;
 const medial=norm(sub(lower,upper));
 const lateralHint=norm([1,0,0]);
 let restNormal=norm(cross(medial,lateralHint));
 if(len(restNormal)<.5)restNormal=[0,0,1];
 const restTangency={};
 for(const name of ['upper','middle','lower']){
  const n=profileNormal(centroid[name]);restTangency[name]=n?Math.abs(dot(restNormal,n)):null;
 }
 return {groups,centroid,restNormal,restTangency,lo,hi};
}
function metrics(model,qS,AC,gapFn,normalFn,scale){
 const regions={},all=[];let globalW=0,globalSq=0,maxOut=-Infinity,maxIn=Infinity;
 for(const name of ['upper','middle','lower']){
  const arr=model.groups[name];let w=0,sq=0,out=-Infinity,inside=Infinity,n=0;
  for(const v of arr){
   if(v.sliceRest==null)continue;
   const p=transform(qS,AC,v.p),g=gapFn(p);if(g==null)continue;
   const d=(g-v.sliceRest)/scale,k=v.weight||1;w+=k;sq+=k*d*d;out=Math.max(out,d);inside=Math.min(inside,d);n++;all.push({d,k});
  }
  const c=transform(qS,AC,sub(model.centroid[name],AC)); // only used if AC is rest-like; replaced below by local centroid calc
  // exact transformed centroid using rest local coordinates:
  const local=sub(model.centroid[name],model.restAC||AC),world=transform(qS,AC,local),tn=normalFn(world);
  const sn=qrot(qS,model.restNormal),tangency=tn?Math.abs(dot(norm(sn),norm(tn))):null,rest=model.restTangency[name];
  regions[name]={n,rms:w?Math.sqrt(sq/w):Infinity,maxOut:out,maxIn:inside,tangency,tangencyDrop:tangency!=null&&rest!=null?Math.max(0,rest-tangency):null};
  globalW+=w;globalSq+=sq;maxOut=Math.max(maxOut,out);maxIn=Math.min(maxIn,inside);
 }
 return {regions,rms:globalW?Math.sqrt(globalSq/globalW):Infinity,maxOut,maxIn};
}
function score(m){
 if(!Number.isFinite(m.rms))return Infinity;
 const cfg={
  upper:{w:1.0,out:10,inward:5,tan:1.0},
  middle:{w:1.35,out:8,inward:4,tan:1.4},
  lower:{w:2.15,out:6,inward:4,tan:2.1}
 };
 let c=.55*m.rms*m.rms;
 for(const name of ['upper','middle','lower']){
  const r=m.regions[name],g=cfg[name];if(!r||!Number.isFinite(r.rms))return Infinity;
  c+=g.w*r.rms*r.rms;
  c+=g.w*3.2*Math.max(0,r.maxOut-g.out)**2;
  c+=g.w*4.3*Math.max(0,-r.maxIn-g.inward)**2;
  if(Number.isFinite(r.tangencyDrop))c+=g.tan*900*r.tangencyDrop*r.tangencyDrop;
 }
 return c;
}
// The old 8/4/2-degree coordinate grid could switch distant clavicle solutions
// between adjacent requested poses. Solve the same weighted residual objective
// jointly, in bounded angle coordinates, without previous-frame state.
function residuals(m){
 const r=[Math.sqrt(.55)*m.rms];
 for(const [name,w,out,inward,tan] of [['upper',1,10,5,1],['middle',1.35,8,4,1.4],['lower',2.15,6,4,2.1]]){
  const g=m.regions[name];
  r.push(Math.sqrt(w)*g.rms,Math.sqrt(w*3.2)*Math.max(0,g.maxOut-out),Math.sqrt(w*4.3)*Math.max(0,-g.maxIn-inward),Math.sqrt(tan*900)*(g.tangencyDrop||0));
 }
 return r;
}
function leastSquares(evaluate,limits,regularizer,iterations=24){
 const n=limits.length;
 let a=limits.map(()=>0),best=evaluate(a),damping=.5,evaluations=1;
 for(let iteration=0;iteration<iterations;iteration++){
  const r=best.residuals;
  if(!r?.every(Number.isFinite))break;
  const J=[],C=best.constraints||[],CJ=[];
  for(let j=0;j<n;j++){
   const lo=a.slice(),hi=a.slice();lo[j]=Math.max(-limits[j],a[j]-.05);hi[j]=Math.min(limits[j],a[j]+.05);
   if(hi[j]===lo[j]){J.push(r.map(()=>0));CJ.push(C.map(()=>0));continue}
   const x=evaluate(lo),y=evaluate(hi);evaluations+=2;
   J.push(r.map((_,k)=>(y.residuals[k]-x.residuals[k])/(hi[j]-lo[j])));
   CJ.push(C.map((_,k)=>(y.constraints[k]-x.constraints[k])/(hi[j]-lo[j])));
  }
  const H=Array.from({length:n},(_,j)=>Array.from({length:2*n+1},(_,k)=>k>n?(k-n-1===j?1:0):k===n?-J[j].reduce((s,v,i)=>s+v*r[i],0)-regularizer[j]*a[j]:J[j].reduce((s,v,i)=>s+v*J[k][i],0)+(j===k?regularizer[j]+damping:0)));
  for(let j=0;j<n;j++){
   let pivot=j;for(let k=j+1;k<n;k++)if(Math.abs(H[k][j])>Math.abs(H[pivot][j]))pivot=k;
   [H[j],H[pivot]]=[H[pivot],H[j]];
   const d=H[j][j];if(Math.abs(d)<1e-12)break;
   for(let k=j;k<=2*n;k++)H[j][k]/=d;
   for(let i=0;i<n;i++)if(i!==j){const v=H[i][j];for(let k=j;k<=2*n;k++)H[i][k]-=v*H[j][k]}
  }
  const delta=H.map(row=>row[n]),planes=C.map((c,k)=>({v:CJ.map(row=>row[k]),b:-c}));
  for(let j=0;j<n;j++)for(const sign of [-1,1]){const v=limits.map(()=>0);v[j]=sign;planes.push({v,b:limits[j]-sign*a[j]})}
  // H-metric projection of the quadratic step onto linearized inequalities.
  // Dual coordinate ascent handles active contacts instead of repeatedly
  // rejecting every step which points across their boundary.
  for(const p of planes){p.h=H.map(row=>p.v.reduce((s,v,j)=>s+row[n+1+j]*v,0));p.d=dot(p.v,p.h);p.lambda=0}
  for(let pass=0;pass<80;pass++){
   let max=0;
   for(const p of planes){if(p.d<1e-12)continue;const error=dot(p.v,delta)-p.b,change=Math.max(0,p.lambda+error/p.d)-p.lambda;p.lambda+=change;for(let j=0;j<n;j++)delta[j]-=p.h[j]*change;max=Math.max(max,Math.abs(change)*Math.sqrt(p.d))}
   if(max<1e-7)break;
  }
  const step=Math.min(1,6/Math.max(Math.hypot(...delta),1e-12));
  let accepted=false;
  for(let line=0;line<8;line++){
   const next=a.map((x,j)=>clamp(x+delta[j]*step*2**(-line),-limits[j],limits[j])),trial=evaluate(next);evaluations++;
   if(trial.total<best.total){a=next;best=trial;accepted=true;break}
  }
  damping=accepted?Math.max(.05,damping*.5):Math.min(1e8,damping*4);
 }
 return {a,best,evaluations};
}
function solve(opts){
 const {baseQC,baseQS,SC,cv,restAC,elev=0,plane=0,points,gapFn,normalFn,scale=1}=opts;
 const identity={qC:baseQC,qS:baseQS,AC:add(SC,qrot(baseQC,cv)),before:null,after:null,params:[0,0,0,0,0,0],metrics:null};
 if(!points?.length||typeof gapFn!=='function'||typeof normalFn!=='function')return identity;
 const model=build(points,restAC,normalFn);if(!model)return identity;model.restAC=restAC.slice();
 const baseAC=add(SC,qrot(baseQC,cv)),baseM=metrics(model,baseQS,baseAC,gapFn,normalFn,scale),baseScore=score(baseM);
 if(!Number.isFinite(baseScore))return {...identity,metrics:baseM};
 const k=smooth(clamp(elev/45,0,1)),limits=[24*k,24*k,20*k,7*k,9*k,8*k];
 const make=a=>{
  const qC0=qm(qa([0,0,1],rad(a[3])),qm(qa([0,1,0],rad(a[4])),baseQC));
  const axis=norm(qrot(qC0,cv)),qC=qm(qa(axis,rad(a[5])),qC0),AC=add(SC,qrot(qC,cv));
  const qS=qm(qa([0,0,1],rad(a[2])),qm(qa([0,1,0],rad(a[0])),qm(qa([1,0,0],rad(a[1])),baseQS)));
  return {qC,qS,AC};
 };
 const evaluate=a=>{
  const p=make(a),m=metrics(model,p.qS,p.AC,gapFn,normalFn,scale),raw=score(m);
  // Corrections remain small and continuous; geometry dominates. Plane-dependent
  // regularizer is deliberately weak because normal kinematics vary substantially.
  const reg=.10*(a[0]**2+a[1]**2)+.08*a[2]**2+.75*(a[3]**2+a[4]**2)+.55*a[5]**2+.002*Math.abs(plane)*a[0]*a[0];
  const constraints=opts.constraints?opts.constraints(p.qS):[];
  return {p,m,raw,total:constraints.some(x=>x>0)?Infinity:raw+reg,residuals:residuals(m),constraints};
 };
 const {a,best,evaluations}=leastSquares(evaluate,limits,[.10+.002*Math.abs(plane),.10,.08,.75,.75,.55]);
 // Never adopt a solution that worsens the geometry-only objective.
 if(best.raw>baseScore+1e-8)return {...identity,before:baseScore,after:baseScore,metrics:baseM};
 return {...best.p,before:baseScore,after:best.raw,params:a,metrics:best.m,model,solver:{method:'bounded-least-squares',evaluations}};
}
function humeralAudit(qS,qH,restShaft){
 const shaft=norm(qrot(qH,restShaft)),local=norm(qrot([-qS[0],-qS[1],-qS[2],qS[3]],shaft));
 // Descriptive only: angle to scapular-plane normal proxy.
 return {shaftWorld:shaft,shaftInScapula:local,planeNormalAngleDeg:Math.acos(clamp(Math.abs(local[2]),-1,1))*180/Math.PI};
}
root.ScapulothoracicEngineV066=Object.freeze({version:'0.66.0',build,metrics,score,solve,humeralAudit,regionOf,residuals,leastSquares});
})(globalThis);
