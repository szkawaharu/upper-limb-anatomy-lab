import '../joint-drag-v071.js';

const E=globalThis.JointDragEngineV071;
const assert=(x,m)=>{if(!x)throw Error(m)};
const rad=d=>d*Math.PI/180;

// Synthetic upper-arm endpoint: elevation + plane describe a spherical limb.
// This is intentionally simple; the test is about gesture mapping stability.
function point(p){
 const e=rad(p.elev),q=rad(p.plane);
 return [Math.sin(e)*Math.cos(q),Math.sin(e)*Math.sin(q),Math.cos(e)];
}
function project(v){return [v[0]*260,-v[1]*260]}
const ranges={elev:[0,180],plane:[0,90],elbow:[0,150]};

// 1) Near the pole: first swipe selects plane exactly once.
const poleAnchor={elev:0,plane:0,axial:23,elbow:0};
const poleLock=E.lockDirection({
 anchor:poleAnchor,current:poleAnchor,side:'right',field:'elev',
 dx:32,dy:-16,point,project,ranges
});
assert(poleLock,'pole lock created');
assert(poleLock.mode==='pole-plane-lock','pole uses one-time plane selection');
const fixedPlane=poleLock.anchor.plane;
let current={...poleAnchor};
for(const [dx,dy] of [[32,-16],[48,-22],[66,-31],[84,-38],[102,-45],[120,-53],[140,-58]]){
 current=E.solve({lock:poleLock,current,dx,dy,point,project,ranges});
 assert(current.plane===fixedPlane,'pole plane must not drift after the initial swipe');
 assert(current.axial===23,'shoulder axial rotation must not be changed');
}
assert(current.elev>20,'pole drag raises the arm');

// 2) Mid-range: lock a single ray in (elevation, plane) space.
// Later pointer wobble may change target position, but it must NOT change the ray.
const midAnchor={elev:58,plane:24,axial:-17,elbow:0};
const midLock=E.lockDirection({
 anchor:midAnchor,current:midAnchor,side:'right',field:'elev',
 dx:18,dy:-11,point,project,ranges
});
assert(midLock&&midLock.mode==='ray','mid-range uses fixed joint-space ray');
const dir=midLock.angleDir.slice();
let mid={...midAnchor},maxRayError=0;
for(const [dx,dy] of [[18,-11],[30,-17],[44,-22],[58,-26],[72,-29],[86,-31],[98,-33]]){
 mid=E.solve({lock:midLock,current:mid,dx,dy,point,project,ranges});
 assert(mid.axial===-17,'axial remains untouched in mid-range drag');
 const de=mid.elev-midLock.anchor.elev,dp=mid.plane-midLock.anchor.plane;
 if(Math.hypot(de,dp)>1){
  // Distance to the locked angular ray; integer sliders permit ~1 degree error.
  const cross=Math.abs(de*dir[1]-dp*dir[0]);
  const denom=Math.max(Math.hypot(...dir),1e-9);
  maxRayError=Math.max(maxRayError,cross/denom);
 }
}
assert(maxRayError<=1.5,'joint-space direction must stay locked');

// 3) Releasing/re-grabbing may create a new mapping (new gesture, new intent).
const second=E.lockDirection({
 anchor:mid,current:mid,side:'right',field:'elev',
 dx:-18,dy:15,point,project,ranges
});
assert(second,'second gesture gets a fresh lock');
assert(
 Math.abs(second.angleDir[0]-midLock.angleDir[0])>.05 ||
 Math.abs(second.angleDir[1]-midLock.angleDir[1])>.05,
 'new gesture may choose a new direction'
);

console.log({
 version:E.version,
 polePlane:fixedPlane,
 poleFinal:{elev:current.elev,plane:current.plane,axial:current.axial},
 midDirection:dir,
 midFinal:{elev:mid.elev,plane:mid.plane,axial:mid.axial},
 maxRayError
});
