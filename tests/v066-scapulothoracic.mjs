import '../scapulothoracic-v066.js';

const E=globalThis.ScapulothoracicEngineV066;
const assert=(x,m)=>{if(!x)throw Error(m)};
const q=[0,0,0,1],SC=[0,0,0],cv=[.45,.82,.46],restAC=cv.slice();

const profileNormal=p=>{
 const v=[p[0],0,p[2]],n=Math.hypot(v[0],v[2])||1;
 return [v[0]/n,0,v[2]/n];
};
const gapFn=p=>Math.hypot(p[0],p[2])-.6;

const points=[];
for(let iy=0;iy<9;iy++)for(let ix=0;ix<7;ix++){
 const y=.5+iy*.08,x=.35+ix*.025,z=Math.sqrt(Math.max(.01,.6*.6-x*x))+.025;
 points.push({p:[x-restAC[0],y-restAC[1],z-restAC[2]],sliceRest:gapFn([x,y,z]),weight:1});
}

let cases=0,nonWorse=0,maxClavicleError=0;
for(const elev of [0,30,60,90,120,150,180])for(const plane of [-90,-45,0,45,90])for(const axial of [-90,0,90]){
 const r=E.solve({baseQC:q,baseQS:q,SC,cv,restAC,elev,plane,points,gapFn,normalFn:profileNormal,scale:.001});
 assert([...r.qC,...r.qS,...r.AC].every(Number.isFinite),'finite');
 const clavicleError=Math.abs(Math.hypot(...r.AC)-Math.hypot(...cv));
 maxClavicleError=Math.max(maxClavicleError,clavicleError);
 assert(clavicleError<1e-9,'SC-AC length');
 if(r.before!==null){assert(r.after<=r.before+1e-8,'geometry non-worsening');nonWorse++}
 if(r.metrics)for(const key of ['upper','middle','lower'])assert(r.metrics.regions[key].n>=2,'three-region samples');
 cases++;
}

const shaft=[0,-1,0],audit=E.humeralAudit(q,q,shaft);
assert(Number.isFinite(audit.planeNormalAngleDeg),'GH descriptive audit');
console.log({version:E.version,cases,nonWorse,maxClavicleError,regions:3,clinicalValidation:false});
