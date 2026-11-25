/* JointDragEngineV0681
 * Stable screen-space inverse kinematics for direct joint dragging.
 * The drag target is anchored to pointer-down pose/screen position so the
 * current articulated pose never becomes a moving target for the next frame.
 */
(function(root){
'use strict';
const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
const add2=(a,b)=>[a[0]+b[0],a[1]+b[1]];
const sub2=(a,b)=>[a[0]-b[0],a[1]-b[1]];
const dot2=(a,b)=>a[0]*b[0]+a[1]*b[1];
const len2=a=>Math.hypot(a[0],a[1]);
const norm2=a=>{const n=Math.max(len2(a),1e-9);return [a[0]/n,a[1]/n]};
const copy=p=>({...p});

function fieldsFor(field){return field==='elbow'?['elbow']:field==='plane'?['plane']:['elev','plane']}
function bounded(v,range){return clamp(v,range[0],range[1])}
function projectPose(pose,field,side,point,project){return project(point(pose,field),side)}

function choosePolePlane(pose,field,side,target,point,project,ranges){
 if(field!=='elev'||!(pose.elev<=2.5||pose.elev>=177.5))return pose;
 const currentXY=projectPose(pose,field,side,point,project),err=sub2(target,currentXY);
 if(len2(err)<3)return pose;
 const dir=norm2(err),[lo,hi]=ranges.plane||[0,90];
 const probeElev=pose.elev<=2.5?Math.min(ranges.elev[1],pose.elev+4):Math.max(ranges.elev[0],pose.elev-4);
 let best=pose.plane,bestScore=-Infinity;
 for(let p=lo;p<=hi+1e-9;p+=5){
  const trial={...pose,elev:probeElev,plane:p},xy=projectPose(trial,field,side,point,project),v=sub2(xy,currentXY),m=len2(v);
  if(m<.05)continue;
  const score=dot2(norm2(v),dir)-Math.abs(p-pose.plane)*0.0008;
  if(score>bestScore){bestScore=score;best=p}
 }
 return {...pose,plane:best};
}

function derivative(pose,f,field,side,point,project,ranges){
 const [lo,hi]=ranges[f],eps=1;
 const a=bounded(pose[f]-eps,[lo,hi]),b=bounded(pose[f]+eps,[lo,hi]);
 if(Math.abs(b-a)<1e-7)return [0,0];
 const pa={...pose,[f]:a},pb={...pose,[f]:b};
 const A=projectPose(pa,field,side,point,project),B=projectPose(pb,field,side,point,project),d=b-a;
 return [(B[0]-A[0])/d,(B[1]-A[1])/d];
}

function targetCost(pose,field,side,target,point,project){
 const q=projectPose(pose,field,side,point,project),d=sub2(q,target);return dot2(d,d);
}

function solve(opts){
 const {anchor,current,side='right',field='elev',dx=0,dy=0,point,project,ranges}=opts||{};
 if(!anchor||!current||typeof point!=='function'||typeof project!=='function'||!ranges)return current;
 const fields=fieldsFor(field),base=projectPose(anchor,field,side,point,project),target=add2(base,[dx,dy]);
 let pose=choosePolePlane(copy(current),field,side,target,point,project,ranges);
 const before=targetCost(current,field,side,target,point,project);
 const xy=projectPose(pose,field,side,point,project),err=sub2(target,xy);
 if(len2(err)<.45)return current;

 const J=fields.map(f=>derivative(pose,f,field,side,point,project,ranges));
 let delta=[];
 if(fields.length===1){
  const j=J[0],den=dot2(j,j)+4;
  delta=[den>1e-8?dot2(j,err)/den:0];
 }else{
  const j0=J[0],j1=J[1];
  const a=dot2(j0,j0)+4,b=dot2(j0,j1),c=dot2(j1,j1)+12;
  const r0=dot2(j0,err),r1=dot2(j1,err),det=a*c-b*b;
  if(Math.abs(det)>1e-8)delta=[(c*r0-b*r1)/det,(a*r1-b*r0)/det];
  else delta=[r0/Math.max(a,1e-8),0];
 }

 const cap=f=>f==='plane'?8:field==='elbow'?7:6;
 let candidate={...pose};
 for(let i=0;i<fields.length;i++){
  const f=fields[i],d=clamp(Number.isFinite(delta[i])?delta[i]:0,-cap(f),cap(f));
  candidate[f]=bounded(candidate[f]+d,ranges[f]);
 }
 for(const f of fields){
  const center=current[f],rounded=Math.round(candidate[f]);
  candidate[f]=bounded(clamp(rounded,center-cap(f),center+cap(f)),ranges[f]);
 }
 if(field==='elev'&&(current.elev<=2.5||current.elev>=177.5))candidate.plane=Math.round(pose.plane);

 let after=targetCost(candidate,field,side,target,point,project);
 if(after>before-.05){
  const half={...current};
  for(const f of fields)half[f]=bounded(Math.round((current[f]+candidate[f])*.5),ranges[f]);
  if(field==='elev'&&(current.elev<=2.5||current.elev>=177.5))half.plane=Math.round(pose.plane);
  const halfCost=targetCost(half,field,side,target,point,project);
  if(halfCost<before-.05){candidate=half;after=halfCost}else return current;
 }
 return candidate;
}
root.JointDragEngineV0681=Object.freeze({version:'0.68.1',solve});
})(globalThis);
