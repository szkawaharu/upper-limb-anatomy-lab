import '../scapulothoracic-v066.js';
import '../scapulohumeral-v067.js';

const V66=globalThis.ScapulothoracicEngineV066,V67=globalThis.ScapulohumeralGuardV067;
const assert=(x,m)=>{if(!x)throw Error(m)};
const norm=a=>{const n=Math.hypot(...a)||1;return a.map(x=>x/n)};
const qa=(axis,a)=>{axis=norm(axis);const s=Math.sin(a/2);return [axis[0]*s,axis[1]*s,axis[2]*s,Math.cos(a/2)]};
const qm=(a,b)=>[a[3]*b[0]+a[0]*b[3]+a[1]*b[2]-a[2]*b[1],a[3]*b[1]-a[0]*b[2]+a[1]*b[3]+a[2]*b[0],a[3]*b[2]+a[0]*b[1]-a[1]*b[0]+a[2]*b[3],a[3]*b[3]-a[0]*b[0]-a[1]*b[1]-a[2]*b[2]];
const qrot=(q,v)=>{const u=q.slice(0,3),cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],add=(a,b)=>a.map((x,i)=>x+b[i]),mul=(a,s)=>a.map(x=>x*s),t=mul(cross(u,v),2);return add(v,add(mul(t,q[3]),cross(u,t)))};
const add=(a,b)=>a.map((x,i)=>x+b[i]);
const SC=[0,0,0],cv=[.45,.82,.46],restAC=cv.slice(),scale=.001;
const gap=p=>Math.hypot(p[0],p[2])-.6;
const normal=p=>{const n=Math.hypot(p[0],p[2])||1;return [p[0]/n,0,p[2]/n]};
const points=[];
for(let iy=0;iy<9;iy++)for(let ix=0;ix<7;ix++){
 const y=.5+iy*.08,x=.35+ix*.025,z=Math.sqrt(Math.max(.01,.36-x*x))+.025;
 points.push({p:[x-restAC[0],y-restAC[1],z-restAC[2]],sliceRest:gap([x,y,z]),weight:1});
}
let cases=0,guarded=0,maxCorrectionDeg=0;
for(const elev of [0,30,60,90,120,150,180])for(const plane of [-90,-45,0,45,90])for(const axial of [-90,0,90]){
 const k=elev/180,f=Math.sin(plane*Math.PI/180);
 const qC=qm(qm(qa([0,0,1],14*k*Math.PI/180),qa([0,1,0],-(18-10*f)*k*Math.PI/180)),qa(cv,25*k*Math.PI/180));
 const AC=add(SC,qrot(qC,cv));
 const qS=qm(qa([0,1,0],(-8+12*f)*k*Math.PI/180),qm(qa([0,0,1],elev/3*Math.PI/180),qa([1,0,0],13*k*Math.PI/180)));
 const phi=plane*Math.PI/180,qH=qm(qa([Math.sin(phi),0,Math.cos(phi)],elev*Math.PI/180),qa([0,-1,0],axial*Math.PI/180));
 const base={qC,qS,AC};
 const baselineAudit=V67.audit(qS,qH,[0,-1,0],elev,plane,qS);
 assert(baselineAudit.proxyOK&&baselineAudit.ghOK&&baselineAudit.correctionOK,'baseline guard');
 const candidate=V66.solve({baseQC:qC,baseQS:qS,SC,cv,restAC,elev,plane,points,gapFn:gap,normalFn:normal,scale});
 const out=V67.refine({v66:V66,base,candidate,SC,cv,restAC,elev,plane,qH,restShaft:[0,-1,0],points,gapFn:gap,normalFn:normal,scale});
 assert([...out.qC,...out.qS,...out.AC].every(Number.isFinite),'finite');
 assert(out.after<=out.before+1e-8,'geometry non-worsening');
 assert(out.v67.proxyOK&&out.v67.ghOK&&out.v67.correctionOK,'guarded result');
 if(out.v67.blend<.999)guarded++;
 maxCorrectionDeg=Math.max(maxCorrectionDeg,out.v67.correctionDeg);
 cases++;
}
console.log({version:V67.version,cases,guarded,maxCorrectionDeg,baselineFalsePositives:0,clinicalValidation:false});
