/* One linked SC->AC->GH model with a simultaneous thoracic-surface fit.
 * No v0.66/0.67/0.72 post-corrections run on this path. All priors are data;
 * all muscles receive the final same rigid transforms. Passive geometry only. */
(function(root){
'use strict';
const M=root.AnatomyMath,F=root.JointFrames,D=root.ShoulderRigData;
const {add,sub,mul,dot,norm,qrot,qmul,qaxis,apply,clamp}=M,rad=x=>x*Math.PI/180;
const thoraxFrame=F.basis([0,0,-1],[0,1,0],[1,0,0]);
// Shape-preserving cubic interpolation with zero endpoint derivative. It
// remains C1 when the educational high-elevation plateau is reached.
const curves=new WeakMap(),surfaces=new WeakMap();
function interpolate(rows,x,col){
 let C=curves.get(rows);
 if(!C){
  C=Array.from({length:rows.length-1},()=>[]);
  for(let c=1;c<rows[0].length;c++){
   const d=rows.slice(0,-1).map((r,i)=>(rows[i+1][c]-r[c])/(rows[i+1][0]-r[0]));
   const m=rows.map((_,i)=>i===0||i===rows.length-1||d[i-1]*d[i]<=0?0:2*d[i-1]*d[i]/(d[i-1]+d[i]));
   for(let i=0;i<rows.length-1;i++){const h=rows[i+1][0]-rows[i][0],a=rows[i][c],b=rows[i+1][c],u=h*m[i],v=h*m[i+1];C[i][c]=[a,u,3*(b-a)-2*u-v,2*(a-b)+u+v]}
  }
  curves.set(rows,C);
 }
 x=clamp(x,rows[0][0],rows.at(-1)[0]);let i=0;while(i<rows.length-2&&x>rows[i+1][0])i++;
 const t=(x-rows[i][0])/(rows[i+1][0]-rows[i][0]),a=C[i][col];return a[0]+t*(a[1]+t*(a[2]+t*a[3]));
}
function prior(elev,plane){
 const p=clamp(plane,0,90),a=p<=40?'abduction':'scapular',b=p<=40?'scapular':'flexion',u=p<=40?p/40:(p-40)/50,w=M.smooth(u);
 return Array.from({length:6},(_,j)=>{const A=D.curves[a],B=D.curves[b];return (interpolate(A,elev,j+1)-A[0][j+1])*(1-w)+(interpolate(B,elev,j+1)-B[0][j+1])*w});
}
function clearance(p,e){
 if(e.rows){
  const rows=e.rows;let table=surfaces.get(e);if(!table){table=rows.map(r=>[r.y,...r.c,...r.r]);surfaces.set(e,table)}
  const cx=interpolate(table,p[1],1),cz=interpolate(table,p[1],2),rx=interpolate(table,p[1],3),rz=interpolate(table,p[1],4),x=p[0]-cx,z=p[2]-cz,r=Math.hypot(x/rx,z/rz);
  return r>1e-9?Math.hypot(x,z)*(1-1/r):-Math.min(rx,rz);
 }
 const v=sub(p,e.c),l=Math.hypot(...v.map((x,j)=>x/e.r[j]));return l>1e-9?Math.hypot(...v)*(1-1/l):-Math.min(...e.r);
}
function linear(A,b){const n=b.length,a=A.map((r,i)=>[...r,b[i]]);for(let j=0;j<n;j++){let k=j;for(let i=j+1;i<n;i++)if(Math.abs(a[i][j])>Math.abs(a[k][j]))k=i;[a[j],a[k]]=[a[k],a[j]];const d=a[j][j];if(!Number.isFinite(d)||Math.abs(d)<1e-12)return null;for(let k=j;k<=n;k++)a[j][k]/=d;for(let i=0;i<n;i++)if(i!==j){const d=a[i][j];for(let k=j;k<=n;k++)a[i][k]-=d*a[j][k]}}return a.map(r=>r[n])}
class ShoulderModel{
 constructor({landmarks,scale=.00555,contact=null,solverIterations=48}){
  this.landmarks=landmarks;this.scale=scale;this.contact=contact;this.solverIterations=solverIterations;
  if(!Number.isFinite(scale)||scale<=0||!Number.isInteger(solverIterations)||solverIterations<1)throw Error('Invalid rig numerical settings');
  const L=landmarks;this.thoraxFrame=L.IJ&&L.C7&&L.PX&&L.T8?F.thorax(L.IJ,L.C7,L.PX,L.T8):thoraxFrame;this.clavicleFrame=F.clavicle(L.SC,L.AC,qrot(this.thoraxFrame,[0,1,0]));this.scapulaFrame=F.scapula(L.TS,L.AA,L.AI);
  this.SC0=F.yxz(qmul(F.conjugate(this.thoraxFrame),this.clavicleFrame));this.AC0=F.yxz(qmul(F.conjugate(this.clavicleFrame),this.scapulaFrame));
  this.neutralH=F.swing(norm(sub(L.EL,L.GH)),qrot(this.thoraxFrame,[0,-1,0]));
  this.cache=new Map();this.girdleCache=new Map();
  this.samples=(contact?.points||[]).map(v=>{const p=v.p||v;return {local:sub(p,L.AC),rest:clearance(p,contact.surface||contact.ellipsoid)}});
  // The proxy is calibrated against the source separation. This preserves the
  // intervening soft-tissue space, rather than forcing bone surfaces together.
  this.sigmaMm=contact?.uncertaintyMm??null;
  this.body=contact&&root.ThoraxContactDataV076&&root.ConvexContact?root.ConvexContact.fromAtlas(root.ThoraxContactDataV076,L.SC,D.landmarks.SC,scale):null;
  this.bodySamples=(this.body?.points||[]).map(p=>sub(p,L.AC));
  this.bodyPlaneCache=new Map();
  this.fieldSignature=root.ShoulderPoseField?.signature({method:'linked-SC-AC-source-tangent-envelope-v4',landmarks,contact,body:this.body,scale,curves:D.curves,prior:D.priorStdDeg,iterations:solverIterations});
  this.field=root.ShoulderPoseFieldData?.signature===this.fieldSignature?root.ShoulderPoseFieldData:null;
  this.fieldReason=this.field?'validated-field':!contact?'contact-data-unavailable':root.ShoulderPoseFieldData?'calibration-signature-mismatch':'field-not-provided';
  if(this.samples.length&&(!Number.isFinite(this.sigmaMm)||this.sigmaMm<=0))throw Error('Contact uncertainty scale must be calibrated');
 }
 forward(params){
  const L=this.landmarks,sc=[this.SC0[0]-rad(params[0]),this.SC0[1]-rad(params[1]),this.SC0[2]+rad(params[2])],ac=[this.AC0[0]+rad(params[3]),this.AC0[1]-rad(params[4]),this.AC0[2]+rad(params[5])];
  const C=qmul(this.thoraxFrame,F.fromYXZ(sc)),S=qmul(C,F.fromYXZ(ac)),qC=qmul(C,F.conjugate(this.clavicleFrame)),qS=qmul(S,F.conjugate(this.scapulaFrame));
  const AC=add(L.SC,qrot(qC,sub(L.AC,L.SC))),GH=add(AC,qrot(qS,sub(L.GH,L.AC)));
  return {qC,qS,AC,GH,SC:L.SC};
 }
 bodyPlanes(target){
  if(!this.body)return [];
  const key=JSON.stringify(target);if(this.bodyPlaneCache.has(key))return this.bodyPlaneCache.get(key);
  const p=this.forward(target),planes=this.bodySamples.map(v=>{const world=add(p.AC,qrot(p.qS,v)),e=root.ConvexContact.evaluate(world,this.body,this.scale,true);return {n:e.gradient,h:dot(e.gradient,world)/this.scale-e.gap};});
  // The smooth support function is convex. Its affine tangent lies below it:
  // a separating tangent halfspace is a conservative bone safeguard. Normals
  // follow the smooth anatomical target, and stay fixed during the joint fit;
  // they never switch between competing envelope faces during optimization.
  if(this.bodyPlaneCache.size>=256)this.bodyPlaneCache.delete(this.bodyPlaneCache.keys().next().value);
  this.bodyPlaneCache.set(key,planes);return planes;
 }
 residual(params,target,includeBody=true){
  const r=params.map((x,j)=>(x-target[j])/D.priorStdDeg[j]);
  // Samples approximate a spatial mean over occupied blade cells; they are
  // not independent anatomical measurements. Squared residuals must therefore
  // integrate to a mean, not grow with discretization density. Each channel
  // (normal clearance and profile domain) uses 1/sqrt(N) in least squares.
  // No anatomical angle, gap limit, prior scale or muscle coefficient changes.
  if(this.samples.length){const p=this.forward(params),denominator=this.sigmaMm*Math.sqrt(this.samples.length);for(const v of this.samples){const d=(clearance(add(p.AC,qrot(p.qS,v.local)),this.contact.surface||this.contact.ellipsoid)-v.rest)/this.scale;r.push(d/denominator);const surface=this.contact.surface,world=add(p.AC,qrot(p.qS,v.local));if(surface?.rows)r.push(Math.max(0,surface.rows[0].y-world[1],world[1]-surface.rows.at(-1).y)/(this.scale*denominator))}}
  if(includeBody&&this.body){const p=this.forward(params),planes=this.bodyPlanes(target),denominator=this.body.residualScaleMm*Math.sqrt(this.bodySamples.length);for(let i=0;i<this.bodySamples.length;i++){const world=add(p.AC,qrot(p.qS,this.bodySamples[i])),plane=planes[i];r.push(Math.min(0,dot(plane.n,world)/this.scale-plane.h)/denominator);}}
  return r;
 }
 jacobian(params,target){
  const r=this.residual(params,target,false),J=[],derivatives=[];
  for(let j=0;j<6;j++){
   const lo=params.slice(),hi=params.slice();lo[j]-=.05;hi[j]+=.05;
   const x=this.residual(lo,target,false),y=this.residual(hi,target,false);
   J.push(r.map((_,i)=>(y[i]-x[i])/.1));
   if(this.body){const a=this.forward(lo),b=this.forward(hi),columns=[[1,0,0],[0,1,0],[0,0,1]].map(v=>mul(sub(qrot(b.qS,v),qrot(a.qS,v)),10));derivatives.push({columns,p:mul(sub(b.AC,a.AC),10)});}
  }
  if(this.body){
   const p=this.forward(params),planes=this.bodyPlanes(target),denominator=this.body.residualScaleMm*Math.sqrt(this.bodySamples.length);
   for(let i=0;i<this.bodySamples.length;i++){
    const v=this.bodySamples[i],plane=planes[i],world=add(p.AC,qrot(p.qS,v)),e={gap:dot(plane.n,world)/this.scale-plane.h,gradient:plane.n};
    r.push(Math.min(0,e.gap)/denominator);
    for(let j=0;j<6;j++){const d=derivatives[j],c=d.columns,dp=[d.p[0]+c[0][0]*v[0]+c[1][0]*v[1]+c[2][0]*v[2],d.p[1]+c[0][1]*v[0]+c[1][1]*v[1]+c[2][1]*v[2],d.p[2]+c[0][2]*v[0]+c[1][2]*v[1]+c[2][2]*v[2]];J[j].push(e.gap<0?dot(e.gradient,dp)/(this.scale*denominator):0);}
   }
  }
  return {r,J};
 }
 pose({elev=0,plane=0,axial=0,sync=true,fit=true,backend="auto"}={}){
  if(![elev,plane,axial].every(Number.isFinite))throw Error('Non-finite joint input');
  const key=JSON.stringify([elev,plane,axial,sync,fit,backend]);if(this.cache.has(key))return this.cache.get(key);
  const target=sync?prior(elev,plane):[0,0,0,0,0,0],gkey=JSON.stringify([elev,plane,sync,fit,backend]),cached=this.girdleCache.get(gkey);let a=cached?.a||target.slice(),lastStep=cached?.lastStep||0;
  // Fixed, damped Gauss-Newton on smooth individual surface residuals. No
  // pose-history seed, grid candidate switching, max-point residual or gate.
  const useField=fit&&sync&&this.field&&backend!=='reference';
  if(!cached&&useField)a=root.ShoulderPoseField.sample(this.field,elev,plane);
  if(!cached&&!useField&&fit&&sync&&this.samples.length)for(let k=0;k<this.solverIterations;k++){
   const {r,J}=this.jacobian(a,target);
   const H=J.map((v,j)=>J.map((w,k)=>v.reduce((s,x,i)=>s+x*w[i],0)+(j===k?.05:0))),g=J.map(v=>-v.reduce((s,x,i)=>s+x*r[i],0)),d=linear(H,g);
   if(!d||!d.every(Number.isFinite))throw Error('Shoulder fit failed');
   let gain=Math.min(1,2/Math.max(Math.hypot(...d),1e-12));
   const loss=r.reduce((s,x)=>s+x*x,0);let next=a;
   // A full GN step may cross several active bone-support constraints.
   // Backtracking prevents an energy-increasing step from producing an
   // unconverged, posture-dependent branch in the precomputed joint field.
   for(let search=0;search<12;search++){
    const candidate=a.map((x,j)=>x+d[j]*gain),residual=this.residual(candidate,target),cost=residual.reduce((s,x)=>s+x*x,0);
    if(cost<=loss){next=candidate;break;}gain*=.5;
   }
   lastStep=next===a?0:Math.hypot(...d)*gain;a=next;
  }
  if(!cached){if(this.girdleCache.size>=256)this.girdleCache.delete(this.girdleCache.keys().next().value);this.girdleCache.set(gkey,{a,lastStep})}
  const p=this.forward(a),L=this.landmarks,shaft=norm(sub(L.EL,L.GH)),axis=qrot(this.thoraxFrame,[-Math.cos(rad(plane)),0,Math.sin(rad(plane))]),qH=qmul(qaxis(axis,rad(elev)),qmul(this.neutralH,qaxis(shaft,rad(axial)))),qGH=qmul(F.conjugate(p.qS),qH);
  const gaps=this.samples.map(v=>(clearance(add(p.AC,qrot(p.qS,v.local)),this.contact.surface||this.contact.ellipsoid)-v.rest)/this.scale);
  const absolute=this.samples.map(v=>clearance(add(p.AC,qrot(p.qS,v.local)),this.contact.surface||this.contact.ellipsoid)/this.scale);
  const outside=this.contact?.surface?.rows?this.samples.filter(v=>{const y=add(p.AC,qrot(p.qS,v.local))[1],r=this.contact.surface.rows;return y<r[0].y||y>r.at(-1).y}).length:0;
  const clavicle={q:p.qC,p:L.SC,rest:L.SC},scapula=F.attach(clavicle,L.AC,qmul(F.conjugate(p.qC),p.qS)),humerus=F.attach(scapula,L.GH,qGH);
  const result={...p,qH:humerus.q,qGH,bones:{clavicle,scapula,humerus},parameters:a,prior:target,assumedScapularUpDeg:null,diagnostics:{version:D.version,method:'linked-SC-AC-GH / calibrated thoracic surface',backend:useField?'validated-field':'reference',fieldReason:this.fieldReason,interpolationErrorBoundMm:useField?this.field.validation.maxRigidErrorBoundMm:null,referenceIterations:this.solverIterations,iterations:useField?0:this.solverIterations,lastIterationStepDeg:lastStep,converged:useField?true:lastStep<.05,contactSamples:gaps.length,surfaceExtrapolatedPoints:outside,absoluteClearanceMm:absolute.length?[Math.min(...absolute),Math.max(...absolute)]:null,clearanceDeltaMm:gaps.length?[Math.min(...gaps),Math.max(...gaps)]:null,clearanceRmsMm:gaps.length?Math.hypot(...gaps)/Math.sqrt(gaps.length):null,clinicalValidation:false,atlasNeutralHumerusOffsetDeg:2*Math.acos(clamp(Math.abs(this.neutralH[3]),0,1))*180/Math.PI,highElevationUnvalidated:elev>120,planeOutsideData:plane<0},model:'SC-AC-GH linked rigid hierarchy; ST surface constraint, literature-informed passive prior'};
  if(![...p.AC,...p.GH,...qH,...p.qS].every(Number.isFinite))throw Error('Non-finite shoulder pose');
  if(this.cache.size>=256)this.cache.delete(this.cache.keys().next().value);this.cache.set(key,result);return result;
 }
}
root.ShoulderKinematics=Object.freeze({ShoulderModel,interpolate,prior,clearance});
})(globalThis);
