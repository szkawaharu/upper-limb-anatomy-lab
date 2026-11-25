import '../joint-drag-v0681.js';

const E=globalThis.JointDragEngineV0681;
const assert=(x,m)=>{if(!x)throw Error(m)};
const rad=x=>x*Math.PI/180;
const point=(p,field)=>field==='elbow'
 ?[p.elbow/100,0,0]
 :[Math.sin(rad(p.elev))*Math.sin(rad(p.plane)),-Math.cos(rad(p.elev)),Math.sin(rad(p.elev))*Math.cos(rad(p.plane))];
const project=p=>[p[0]*300,p[1]*300];
const ranges={elev:[0,180],plane:[0,90],elbow:[0,150]};

let current={elev:0,plane:0,elbow:0},anchor={...current};
let visibleMaxStep=0,signFlips=0,lastElevStep=0;
for(let i=1;i<=120;i++){
 const next=E.solve({anchor,current,field:'elev',dx:i*1.7,dy:-i*.8,point,project,ranges});
 const de=next.elev-current.elev,dp=next.plane-current.plane;
 if(current.elev>2.5&&current.elev<177.5)visibleMaxStep=Math.max(visibleMaxStep,Math.abs(de),Math.abs(dp));
 if(lastElevStep&&de&&Math.sign(lastElevStep)!==Math.sign(de))signFlips++;
 lastElevStep=de||lastElevStep;
 assert(Number.isFinite(next.elev)&&Number.isFinite(next.plane),'finite shoulder pose');
 current=next;
}
assert(visibleMaxStep<=8,'bounded visible per-event step');
assert(signFlips<=2,'no alternating shoulder feedback oscillation');

// Pointer held still after convergence must not create a two-state cycle.
const heldDx=204,heldDy=-96,seen=[];
for(let i=0;i<30;i++){
 current=E.solve({anchor,current,field:'elev',dx:heldDx,dy:heldDy,point,project,ranges});
 seen.push(current.elev+':'+current.plane);
}
assert(new Set(seen.slice(-10)).size<=2,'held pointer converges instead of trembling');

// Elbow remains a monotone one-DOF solve under monotone drag.
current={elev:70,plane:35,elbow:20};anchor={...current};let last=current.elbow;
for(let i=1;i<=50;i++){
 current=E.solve({anchor,current,field:'elbow',dx:i*2,dy:0,point,project,ranges});
 assert(Math.abs(current.elbow-last)<=7,'elbow per-event cap');last=current.elbow;
}
console.log({version:E.version,visibleMaxStep,signFlips,heldStates:new Set(seen.slice(-10)).size,finalShoulder:seen.at(-1),finalElbow:current.elbow});
