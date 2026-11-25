import '../scapulothoracic-v066.js';
import '../scapulohumeral-v067.js';
import '../scapulothoracic-v072.js';

const V66=globalThis.ScapulothoracicEngineV066;
const V67=globalThis.ScapulohumeralGuardV067;
const V72=globalThis.ScapulothoracicContactGuardV072;
const assert=(x,m)=>{if(!x)throw Error(m)};
const norm=a=>{const n=Math.hypot(...a)||1;return a.map(x=>x/n)};
const qa=(axis,a)=>{axis=norm(axis);const s=Math.sin(a/2);return [axis[0]*s,axis[1]*s,axis[2]*s,Math.cos(a/2)]};
const qm=(a,b)=>[a[3]*b[0]+a[0]*b[3]+a[1]*b[2]-a[2]*b[1],a[3]*b[1]-a[0]*b[2]+a[1]*b[3]+a[2]*b[0],a[3]*b[2]+a[0]*b[1]-a[1]*b[0]+a[2]*b[3],a[3]*b[3]-a[0]*b[0]-a[1]*b[1]-a[2]*b[2]];
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const add=(a,b)=>a.map((x,i)=>x+b[i]),mul=(a,s)=>a.map(x=>x*s);
const qrot=(q,v)=>{const u=q.slice(0,3),t=mul(cross(u,v),2);return add(v,add(mul(t,q[3]),cross(u,t)))};

const SC=[0,0,0],cv=[.45,.82,.46],restAC=cv.slice(),scale=.001;
const gap=p=>Math.hypot(p[0],p[2])-.6;
const normal=p=>{const n=Math.hypot(p[0],p[2])||1;return [p[0]/n,0,p[2]/n]};
const points=[];
for(let iy=0;iy<9;iy++)for(let ix=0;ix<7;ix++){
 const y=.5+iy*.08,x=.35+ix*.025,z=Math.sqrt(Math.max(.01,.36-x*x))+.025;
 points.push({p:[x-restAC[0],y-restAC[1],z-restAC[2]],sliceRest:gap([x,y,z]),weight:1});
}

function run(elev,plane,axial=0){
 const k=elev/180,f=Math.sin(plane*Math.PI/180);
 const qC=qm(qm(qa([0,0,1],14*k*Math.PI/180),qa([0,1,0],-(18-10*f)*k*Math.PI/180)),qa(cv,25*k*Math.PI/180));
 const AC=add(SC,qrot(qC,cv));
 const qS=qm(qa([0,1,0],(-8+12*f)*k*Math.PI/180),qm(qa([0,0,1],elev/3*Math.PI/180),qa([1,0,0],13*k*Math.PI/180)));
 const phi=plane*Math.PI/180,qH=qm(qa([Math.sin(phi),0,Math.cos(phi)],elev*Math.PI/180),qa([0,-1,0],axial*Math.PI/180));
 const base={qC,qS,AC};
 const candidate=V66.solve({baseQC:qC,baseQS:qS,SC,cv,restAC,elev,plane,points,gapFn:gap,normalFn:normal,scale});
 const guarded=V67.refine({v66:V66,base,candidate,SC,cv,restAC,elev,plane,qH,restShaft:[0,-1,0],points,gapFn:gap,normalFn:normal,scale});
 return V72.refine({v66:V66,v67:V67,input:guarded,base,SC,cv,restAC,elev,plane,qH,restShaft:[0,-1,0],points,gapFn:gap,normalFn:normal,scale});
}

let cases=0,active=0,regionalRegressions=0,maxCorrection=0;
for(const elev of [30,60,90,120,150,180])for(const plane of [-90,-45,0,45,90]){
 const o=run(elev,plane);
 assert([...o.qC,...o.qS,...o.AC].every(Number.isFinite),'finite');
 assert(o.v72&&Number.isFinite(o.v72.beforeScore)&&Number.isFinite(o.v72.afterScore),'v72 diagnostics');
 if(o.v72.active){
  active++;
  assert(o.v72.afterScore<o.v72.beforeScore,'active correction improves contact objective');
  maxCorrection=Math.max(maxCorrection,o.v72.correctionDeg);
 }
 for(const k of ['upper','middle','lower']){
  const a=V72.violation(o.v72.after.regions[k],V72.CFG[k]);
  const b=V72.violation(o.v72.before.regions[k],V72.CFG[k]);
  const bad=a.out>b.out+.500001||a.inside>b.inside+.500001||a.rms>b.rms+.500001||a.tan>b.tan+.020001;
  if(bad)regionalRegressions++;
  assert(!bad,'no regional contact tradeoff '+elev+' '+plane+' '+k);
 }
 if(o.v72.audit)assert(V72.hardKinematicOK(o.v72.audit),'expanded kinematic guard');
 cases++;
}
assert(active>0,'guard activates in synthetic winging/penetration cases');
assert(regionalRegressions===0,'no region becomes materially worse');

// Specific case where v0.67 backs away from the v0.66 contact solution.
const example=run(60,0);
assert(example.v72.active,'known contact-loss case activates');
assert(example.v72.after.regions.lower.maxOut<example.v72.before.regions.lower.maxOut-.5,'lower-border separation improves');

console.log({
 version:V72.version,cases,active,regionalRegressions,maxCorrection,
 exampleLowerBefore:example.v72.before.regions.lower.maxOut,
 exampleLowerAfter:example.v72.after.regions.lower.maxOut,
 clinicalValidation:false
});
