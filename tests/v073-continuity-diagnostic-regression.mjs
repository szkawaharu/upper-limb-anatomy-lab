// Known-defect reproducer. A passing test here means the defect was detected,
 // NOT that the production frame is continuous.
import assert from "node:assert/strict";
import {trackStep} from "./browser-step-diagnostics.mjs";
const M=globalThis.AnatomyMath;
const t=[1,0,0],u=[0,1,0],near=[.0142,.0141].map(x=>[-Math.cos(x),Math.sin(x),0]);
const outputs=near.map(v=>M.align(t,v,u,u));
const jump=M.len(M.sub(outputs[0],outputs[1]));
assert(jump>1.99,"old alignment reproducer");
const P=globalThis.FiberPath,restFrame={c:[0,0,0],t,u,n:[0,0,1],tl:1,ul:1},anchor=M.identity();
const frames=near.map(v=>P.transported(restFrame,{...restFrame,t:v},anchor,anchor,.5));
const fixedJump=M.len(M.sub(frames[0].u,frames[1].u));
assert(fixedJump<1e-4,"no hard-switch jump");
assert(frames.every(f=>Math.abs(M.dot(M.cross(f.t,f.u),f.n)-1)<1e-9));
assert(frames.some(f=>f.frameConfidence<1&&f.effectiveTangentDot<0),"ill-conditioned tangent explicitly reported");
for(let i=0;i<=1000;i++){
 const theta=i*Math.PI/1000,v=[Math.cos(theta),Math.sin(theta),0],f=P.transported(restFrame,{...restFrame,t:v},anchor,anchor,.5);
 assert([f.t,f.u,f.n].flat().every(Number.isFinite));
 assert(Math.abs(M.dot(M.cross(f.t,f.u),f.n)-1)<1e-9);
 if(f.frameConfidence===1)assert(M.dot(f.t,v)>1-1e-9);
}
const a=M.identity(),b=M.identity(),rest=[[[0,0,0],[1/3,0,0],[2/3,0,0],[1,0,0]],[[0,.2,0],[1/3,.2,0],[2/3,.2,0],[1,.2,0]]];
const state={a,b,restControls:rest,posedControls:rest,sheet:true,blendRange:[0,1],tendonStart:.8,reference:[0,1,0],contact:null};
const row={id:"synthetic",name:"synthetic",side:"right",state,samples:[{vertexId:42,p:[.5,.1,.03],q:[.5,.1,.03],t:.5,l:.5,materialT:.5}]};
const maxima={};
assert.equal(trackStep(maxima,row,row,0,{plane:0,fromElev:0,toElev:3,direction:"forward"}),0);
assert.equal(maxima.synthetic.vertexId,42);
assert(maxima.synthetic.after.reconstructedErrorMm<1e-9);
assert.equal(maxima.synthetic.after.tangentDot,1);
assert.throws(()=>trackStep(maxima,row,{...row,samples:[{...row.samples[0],vertexId:43}]},0,{}));
console.log(JSON.stringify({knownDefect:"antiparallel frame axis jump",angleStepDeg:.0001*180/Math.PI,jumpUnitVector:jump,candidateJumpUnitVector:fixedJump,diagnosticIdentity:true,qualityPass:false}));
