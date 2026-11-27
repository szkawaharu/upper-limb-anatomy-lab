import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
await import('../engine/math.js');await import('../engine/camera-framing.js');
const C=CameraFraming,box=(lo,hi)=>[lo[0],hi[0]].flatMap(x=>[lo[1],hi[1]].flatMap(y=>[lo[2],hi[2]].map(z=>[x,y,z])));
const neutral=[...box([-.75,-.05,-.36],[.75,1.85,.38]),...box([-.19,1.85,-.25],[.19,2.25,.24]),...box([-1.02,-1.65,-.13],[-.75,1.4,.14]),...box([.75,-1.65,-.13],[1.02,1.4,.14])];
// Synthetic display-bound boxes only, not a new anatomical mesh.
const raised=[...neutral.slice(0,16),...box([-.98,1.3,-.15],[-.67,4.22,.16]),...neutral.slice(24)];
const both=[...neutral.slice(0,16),...raised.slice(16,24),...box([.75,1.25,-.5],[1.72,3.83,.14])];
let tested=0;
for(const [w,h] of [[320,568],[390,844],[512,850],[768,1024],[1024,768]]){
 const pads=[Math.min(44,Math.max(20,w*.07)),Math.min(72,Math.max(52,w*.16)),Math.min(96,Math.max(64,h*.08)),Math.min(92,Math.max(56,h*.075))];
 for(const yaw of [0,Math.PI/2,Math.PI,Math.PI*.75])for(const pitch of [-.35,.12,.5]){
  const initial={yaw,pitch,target:[0,0,0],zoom:4};
  const base=C.fit(neutral,initial,w,h,pads);
  assert(!C.clipped(neutral,base,w,h,pads));
  assert.equal(C.protect(neutral,base,w,h,pads),null);
  for(const points of [raised,both]){
   const next=C.protect(points,base,w,h,pads)||base;
   assert(next.zoom>=base.zoom,'guard must never zoom in');
   assert.equal(next.yaw,yaw);assert.equal(next.pitch,pitch);
   assert(!C.clipped(points,next,w,h,pads),'all posed bounds in usable viewport');
   assert.equal(C.protect(points,next,w,h,pads),null,'idempotent');
   const rotated={...next,yaw:yaw+Math.PI*.63};
   const fit=C.protect(points,rotated,w,h,pads)||rotated;
   assert(!C.clipped(points,fit,w,h,pads));
   const manual={...fit,zoom:fit.zoom*.7,target:[.2,-.2,.4]};
   assert.deepEqual(C.screenBounds(points,manual,w,h),C.screenBounds(points,manual,w,h),'diagnosis is pure');
   tested++;
  }
  // Wider previous span must be used in the centering calculation itself.
  const wide=C.fit(raised,base,w,h,pads,{minSpan:base.zoom*2});
  const s=C.screenBounds(raised,wide,w,h),l=C.layout(w,h,pads);
  assert(Math.abs((s.minX+s.maxX)/2-l.cx)<1e-8);
  assert(Math.abs((s.minY+s.maxY)/2-l.cy)<1e-8);
 }
}
const html=readFileSync(new URL('../index.html',import.meta.url),'utf8');
assert(/id="v572Fit" type="checkbox" checked/.test(html),'clipping guard default enabled');
assert(html.includes('CameraFraming.protect(pts,state,w,h'));
assert(html.includes('if(pointers.size||v73NeckDrag)return;v572FinishInteraction()'),'idle must not refit during active drag');
assert(html.includes('v572ManualFraming=true;v572OrbitNeedsFit=false;'),'manual pinch/pan cancels orbit guard');
assert(html.includes('v572NeedsFit||(v572OrbitNeedsFit&&!v572ManualFraming)'));
console.log(JSON.stringify({tested,viewports:5,views:12,posedCases:2,syntheticDisplayBounds:true,anatomicalAcceptance:false,browserDeviceTest:false}));
