/* JointDragEngineV071
 * Direction-locked shoulder drag for direct manipulation.
 *
 * The first meaningful swipe establishes ONE mapping from screen direction to
 * the shoulder's (elevation, elevation-plane) path. That mapping is held until
 * pointer-up. The current pose is never allowed to redefine what the same
 * swipe direction means midway through the gesture.
 *
 * This is a UI/kinematic constraint, not a biomechanical shoulder model.
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

function bounded(v,range){return clamp(v,range[0],range[1])}
function projectPose(pose,field,side,point,project){return project(point(pose,field),side)}

function derivative(pose,f,field,side,point,project,ranges,eps=1){
 const [lo,hi]=ranges[f],a=bounded(pose[f]-eps,[lo,hi]),b=bounded(pose[f]+eps,[lo,hi]);
 if(Math.abs(b-a)<1e-8)return [0,0];
 const A=projectPose({...pose,[f]:a},field,side,point,project),B=projectPose({...pose,[f]:b},field,side,point,project),d=b-a;
 return [(B[0]-A[0])/d,(B[1]-A[1])/d];
}

function poleLock(anchor,side,dir,point,project,ranges){
 const base=projectPose(anchor,'elev',side,point,project),[elo,ehi]=ranges.elev,[plo,phi]=ranges.plane;
 let best=null;
 const probe=12;
 for(const sign of [-1,1]){
  const e=bounded(anchor.elev+sign*probe,[elo,ehi]);
  if(Math.abs(e-anchor.elev)<2)return best;
  for(let plane=plo;plane<=phi+1e-9;plane+=1){
   const q=projectPose({...anchor,elev:e,plane},'elev',side,point,project),v=sub2(q,base),m=len2(v);
   if(m<.05)continue;
   const score=dot2(norm2(v),dir)-Math.abs(plane-anchor.plane)*0.00002;
   if(!best||score>best.score)best={score,plane,sign};
  }
 }
 return best;
}

function rayFromJacobian(anchor,side,dir,point,project,ranges){
 const jE=derivative(anchor,'elev','elev',side,point,project,ranges,1.5);
 const jP=derivative(anchor,'plane','elev',side,point,project,ranges,1.5);
 // Damped least squares. Plane receives stronger damping because at high/low
 // elevation it can mimic an axial-looking roll in screen space.
 const a=dot2(jE,jE)+2.5,b=dot2(jE,jP),c=dot2(jP,jP)+9.0;
 const r0=dot2(jE,dir),r1=dot2(jP,dir),det=a*c-b*b;
 let de=0,dp=0;
 if(Math.abs(det)>1e-9){de=(c*r0-b*r1)/det;dp=(a*r1-b*r0)/det}
 else de=dot2(jE,dir)/(dot2(jE,jE)+2.5);
 if(!Number.isFinite(de))de=0;if(!Number.isFinite(dp))dp=0;
 const m=Math.max(Math.abs(de),Math.abs(dp));
 if(m<1e-7){
  const s=dot2(jE,dir)>=0?1:-1;
  return [s,0];
 }
 return [de/m,dp/m];
}

function lambdaRange(anchor,angleDir,ranges){
 let lo=-Infinity,hi=Infinity;
 for(const [i,f] of ['elev','plane'].entries()){
  const d=angleDir[i];if(Math.abs(d)<1e-9)continue;
  const a=(ranges[f][0]-anchor[f])/d,b=(ranges[f][1]-anchor[f])/d;
  lo=Math.max(lo,Math.min(a,b));hi=Math.min(hi,Math.max(a,b));
 }
 if(!Number.isFinite(lo))lo=-180;if(!Number.isFinite(hi))hi=180;
 return lo<=hi?[lo,hi]:[0,0];
}

function lockDirection(opts){
 const {anchor,side='right',field='elev',dx=0,dy=0,point,project,ranges,minPixels=7}=opts||{};
 if(field!=='elev'||!anchor||typeof point!=='function'||typeof project!=='function'||!ranges?.elev||!ranges?.plane)return null;
 const mag=Math.hypot(dx,dy);if(!Number.isFinite(mag)||mag<minPixels)return null;
 const dir=norm2([dx,dy]),base=projectPose(anchor,field,side,point,project);
 let lockedAnchor=copy(anchor),angleDir,mode='ray',pole=null;
 // Plane is poorly defined close to the arm-down / arm-up pole. Resolve it
 // ONCE from the first swipe direction, then freeze that plane for this drag.
 if(anchor.elev<=10||anchor.elev>=170){
  pole=poleLock(anchor,side,dir,point,project,ranges);
  if(pole){
   lockedAnchor.plane=pole.plane;
   angleDir=[pole.sign,0];
   mode='pole-plane-lock';
  }
 }
 if(!angleDir)angleDir=rayFromJacobian(lockedAnchor,side,dir,point,project,ranges);
 const range=lambdaRange(lockedAnchor,angleDir,ranges);
 return {
  version:'0.71.0',mode,field,side,
  anchor:lockedAnchor,
  anchorScreen:base.slice(),
  initialScreenDir:dir.slice(),
  angleDir:angleDir.slice(),
  lambda:0,
  lambdaRange:range,
  polePlane:pole?.plane??null
 };
}

function poseAt(lock,lambda,ranges){
 const p={...lock.anchor};
 p.elev=bounded(lock.anchor.elev+lock.angleDir[0]*lambda,ranges.elev);
 p.plane=bounded(lock.anchor.plane+lock.angleDir[1]*lambda,ranges.plane);
 // Sliders are integer degrees; keeping rounding here makes CPU/UI state and
 // rendered state identical and avoids sub-degree feedback oscillation.
 p.elev=Math.round(p.elev);p.plane=Math.round(p.plane);
 return p;
}

function targetCost(pose,side,target,point,project){
 const q=projectPose(pose,'elev',side,point,project),d=sub2(q,target);
 return dot2(d,d);
}

function solve(opts){
 const {lock,current,dx=0,dy=0,point,project,ranges}=opts||{};
 if(!lock||!current||lock.field!=='elev'||typeof point!=='function'||typeof project!=='function'||!ranges)return current;
 const target=add2(lock.anchorScreen,[dx,dy]),[globalLo,globalHi]=lock.lambdaRange;
 let center=clamp(Number.isFinite(lock.lambda)?lock.lambda:0,globalLo,globalHi);
 // Keep one event local. This is the second continuity guard after the fixed
 // gesture ray: sparse touch events cannot jump to a remote projected solution.
 const lo=Math.max(globalLo,center-8),hi=Math.min(globalHi,center+8);
 const currentCost=targetCost(current,lock.side,target,point,project);
 let bestPose=current,bestLambda=center,bestCost=currentCost;
 const tryLambda=L=>{
  const l=clamp(L,lo,hi),p=poseAt(lock,l,ranges);
  const c=targetCost(p,lock.side,target,point,project)+0.012*(l-center)*(l-center);
  if(c<bestCost-.02){bestCost=c;bestPose=p;bestLambda=l}
 };
 for(const d of [-8,-6,-4,-2,0,2,4,6,8])tryLambda(center+d);
 for(const step of [1,.5]){
  const c=bestLambda;for(const d of [-step,step])tryLambda(c+d);
 }
 if(bestPose===current)return current;
 lock.lambda=bestLambda;
 // Never modify axial rotation, pronation, elbow, etc. in a shoulder drag.
 return {...current,elev:bestPose.elev,plane:bestPose.plane};
}

root.JointDragEngineV071=Object.freeze({version:'0.71.0',lockDirection,solve});
})(globalThis);
