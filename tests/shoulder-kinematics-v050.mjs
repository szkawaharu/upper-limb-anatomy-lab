// Usage: node tests/shoulder-kinematics-v050.mjs (repository root)
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';

const html=readFileSync('index.html','utf8');
const js=html.match(/<script>([\s\S]*?)<\/script>/)?.[1];
assert.ok(js, 'inline JS source missing');
new Function(js);
const start=js.indexOf('function shoulderPoseFor(');
const end=js.indexOf('function makeRealPose(){',start);
assert.ok(start>0&&end>start,'pure shoulder evaluator missing');
const add=(a,b)=>a.map((x,i)=>x+b[i]);
const sub=(a,b)=>a.map((x,i)=>x-b[i]);
const mul=(a,k)=>a.map(x=>x*k);
const dot=(a,b)=>a.reduce((v,x,i)=>v+x*b[i],0);
const len=a=>Math.hypot(...a);
const norm=a=>mul(a,1/Math.max(len(a),1e-12));
const clamp=(v,a,b)=>Math.max(a,Math.min(b,v));
const rad=d=>d*Math.PI/180;
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const qa=(axis,r)=>[...mul(norm(axis),Math.sin(r/2)),Math.cos(r/2)];
const qm=(a,b)=>[a[3]*b[0]+a[0]*b[3]+a[1]*b[2]-a[2]*b[1],a[3]*b[1]-a[0]*b[2]+a[1]*b[3]+a[2]*b[0],a[3]*b[2]+a[0]*b[1]-a[1]*b[0]+a[2]*b[3],a[3]*b[3]-a[0]*b[0]-a[1]*b[1]-a[2]*b[2]];
const qc=q=>[-q[0],-q[1],-q[2],q[3]];
const qr=(q,v)=>{const t=mul(cross(q.slice(0,3),v),2);return add(v,add(mul(t,q[3]),cross(q.slice(0,3),t)))};
const transform=(q,p,v)=>add(qr(q,v),p);
const qt=(from,to)=>{const a=norm(from),b=norm(to),d=dot(a,b);
 if(d<-.999999){let axis=cross([1,0,0],a);if(len(axis)<1e-5)axis=cross([0,0,1],a);return qa(axis,Math.PI)}
 const v=[...cross(a,b),1+d],m=Math.hypot(...v);return v.map(x=>x/m)};
const realAnchors={SC:[0,0,0],AC:[.55,.02,.12],GH:[.65,-.16,.15],EL:[.95,-1.6,.25]};
const settings={elev:0,plane:0,axial:0,sync:true},atlasScale=.00555,ID=[0,0,0,1];
const api={add,sub,mul,dot,len,norm,clamp,rad,qa,qm,qc,qr,qt,transform,realAnchors,settings,atlasScale,ID};
const pure=js.slice(start,end);
const {shoulderPoseFor,shoulderPoseMetrics}=new Function(...Object.keys(api),pure+'; return {shoulderPoseFor,shoulderPoseMetrics}')(...Object.values(api));
let count=0;
const baseLen=len(sub(realAnchors.AC,realAnchors.SC))/atlasScale;
for(const sync of [false,true])for(const elev of [0,30,60,90,120,150,180])for(const plane of [0,45,90])for(const axial of [-45,0,45]){
 const m=shoulderPoseMetrics(elev,plane,axial,sync);
 assert.ok([m.scGapMm,m.acGapMm,m.ghGapMm,m.clavicleLengthMm,m.ghDisplacementMm,m.relativeTiltDeg].every(Number.isFinite));
 assert.ok(Math.max(m.scGapMm,m.acGapMm,m.ghGapMm)<1e-8,'shoulder anchor continuity');
 assert.ok(Math.abs(m.clavicleLengthMm-baseLen)<1e-8,'clavicle length invariant');
 assert.ok(m.relativeTiltDeg>=0&&m.relativeTiltDeg<180.01,'relative tilt finite');
 assert.equal(m.clinicalAccuracyVerified,false);
 const p=shoulderPoseFor(elev,plane,axial,sync);
 const k=elev/180,yaw=rad(sync?12*k:0);
 const qOld=qm(qa([0,0,1],rad(sync?elev/3:0)),qm(qa([1,0,0],rad(sync?13*k:0)),qa([0,1,0],rad(sync?10*k:0))));
 const qScap=qm(qa([0,1,0],yaw),qOld);
 const qClav=sync?qm(qa([0,0,1],rad(18*k)),qa([0,1,0],yaw)):ID;
 const theta=rad(elev),phi=rad(plane),dir=elev===180?[.001,1,.001]:[Math.sin(theta)*Math.cos(phi),-Math.cos(theta),Math.sin(theta)*Math.sin(phi)];
 const qHum=qm(qt([0,-1,0],dir),qa([0,-1,0],rad(axial)));
 const eq=(a,b)=>len(sub(a,b))<1e-10;
 assert.ok(eq(qScap,p.qS)&&eq(qClav,p.qC)&&eq(qHum,p.qH),'existing kinematics changed');
 count++;
}
const ids=['shoulderAuditPanel','shoulderAuditCurrent','shoulderAuditSweep','shoulderAuditSweepResult'];
for(const id of ids)assert.ok(html.includes('id="'+id+'"'),'missing UI: '+id);
assert.ok(html.includes('上肢3D解剖ラボ · v0.55'),'version 0.51 not visible');
assert.ok(js.includes('window.__labShoulder=()=>'),'missing user inspect API');
console.log('PASS v0.50: '+count+' posture combinations, SC/AC/GH continuity and legacy motion preservation. Clinical accuracy is NOT certified.');
