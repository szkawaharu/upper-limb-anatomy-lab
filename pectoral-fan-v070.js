/* PectoralFanEngineV070
 * Pectoralis-major-specific 2-D fan-sheet deformation.
 * Six transverse fascicle lanes span the broad origin and converge toward the
 * humeral insertion. A local sheet frame maps the atlas surface into the posed
 * fan. Educational geometry only; no activation/force/contact/FEM.
 */
(function(root){
'use strict';
const N=6,clamp=(x,a,b)=>Math.max(a,Math.min(b,x)),smooth=x=>{x=clamp(x,0,1);return x*x*(3-2*x)};
const add=(a,b)=>a.map((x,i)=>x+b[i]),sub=(a,b)=>a.map((x,i)=>x-b[i]),mul=(a,s)=>a.map(x=>x*s),dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],len=a=>Math.hypot(...a),norm=a=>mul(a,1/Math.max(len(a),1e-12)),mix=(a,b,t)=>add(mul(a,1-t),mul(b,t));
const qrot=(q,v)=>{const u=q.slice(0,3),t=mul(cross(u,v),2);return add(v,add(mul(t,q[3]),cross(u,t)))},qmix=(a,b,t)=>norm(mix(a,mul(b,dot(a,b)<0?-1:1),t)),apply=(a,p)=>add(a.p,qrot(a.q,sub(p,a.rest)));
function curve(c,t){const u=1-t;return add(add(mul(c[0],u*u*u),mul(c[1],3*u*u*t)),add(mul(c[2],3*u*t*t),mul(c[3],t*t*t)))}
function tangent(c,t){const u=1-t;return add(add(mul(sub(c[1],c[0]),3*u*u),mul(sub(c[2],c[1]),6*u*t)),mul(sub(c[3],c[2]),3*t*t))}
function laneCurve(controls,l){const x=clamp(l,0,N-1),i=Math.min(N-2,Math.floor(x)),f=x-i;return controls[i].map((p,j)=>mix(p,controls[i+1][j],f))}
function mean(points,weights){let s=[0,0,0],d=0;for(let i=0;i<points.length;i++){const w=weights[i];s=add(s,mul(points[i],w));d+=w}return d>1e-12?mul(s,1/d):points[0].slice()}
function percentile(a,p){const s=a.slice().sort((x,y)=>x-y);return s[Math.max(0,Math.min(s.length-1,Math.floor(p*(s.length-1))))]}
function transverseAxis(vertices,weights,O,I){
 const d=norm(sub(I,O));let v=sub([0,1,0],mul(d,dot([0,1,0],d)));if(len(v)<.1)v=sub([0,0,1],mul(d,dot([0,0,1],d)));v=norm(v);
 const C=[[0,0,0],[0,0,0],[0,0,0]];
 for(let i=0;i<vertices.length;i++){const p=vertices[i],t=clamp(weights[i],0,1),line=mix(O,I,t);let r=sub(p,line);r=sub(r,mul(d,dot(r,d)));const w=.15+(1-t)*(1-t);for(let a=0;a<3;a++)for(let b=0;b<3;b++)C[a][b]+=r[a]*r[b]*w}
 for(let k=0;k<12;k++){let q=C.map(row=>dot(row,v));q=sub(q,mul(d,dot(q,d)));if(len(q)<1e-12)break;v=norm(q)}
 if(Math.abs(dot(v,[0,1,0]))>.12&&dot(v,[0,1,0])<0)v=mul(v,-1);
 return {d,u:v,n:norm(cross(d,v))};
}
function laneValue(p,t,field){const r=sub(p,mix(field.origin,field.insertion,clamp(t,0,1))),x=dot(r,field.u);return clamp((x-field.uLo)/Math.max(field.uHi-field.uLo,1e-8),0,1)*(N-1)}
function build(vertices,faces,weights){
 if(!vertices?.length||vertices.length!==weights?.length||vertices.length<24||!faces?.length)return null;
 const O=mean(vertices,weights.map(w=>(1-clamp(w,0,1))**4+.001)),I=mean(vertices,weights.map(w=>clamp(w,0,1)**4+.001));if(len(sub(I,O))<.03)return null;
 const axes=transverseAxis(vertices,weights,O,I),raw=vertices.map((p,i)=>dot(sub(p,mix(O,I,clamp(weights[i],0,1))),axes.u)),broad=raw.filter((_,i)=>weights[i]<.82);
 const uLo=percentile(broad.length>20?broad:raw,.025),uHi=percentile(broad.length>20?broad:raw,.975);if(!Number.isFinite(uLo)||!Number.isFinite(uHi)||uHi-uLo<.015)return null;
 const field={engine:'fan070',origin:O,insertion:I,u:axes.u,d:axes.d,n:axes.n,uLo,uHi},uv=new Float32Array(vertices.length);
 const cells=Array.from({length:N},()=>Array.from({length:4},()=>({s:[0,0,0],w:0}))),global=Array.from({length:4},()=>({s:[0,0,0],w:0}));
 const addCell=(c,p,w)=>{if(w>0){c.s=add(c.s,mul(p,w));c.w+=w}};
 for(let i=0;i<vertices.length;i++){
  const p=vertices[i],t=clamp(weights[i],0,1),l=laneValue(p,t,field);uv[i]=l;const st=t*3;
  for(let j=Math.max(0,Math.floor(st)-1);j<=Math.min(3,Math.ceil(st)+1);j++){
   const wt=Math.max(0,1-Math.abs(j-st)/1.15)**2;addCell(global[j],p,wt);
   for(let k=Math.max(0,Math.floor(l)-1);k<=Math.min(N-1,Math.ceil(l)+1);k++)addCell(cells[k][j],p,wt*Math.max(0,1-Math.abs(k-l)/1.25)**2);
  }
 }
 if(global.some(c=>c.w<1e-8))return null;const gc=global.map(c=>mul(c.s,1/c.w));
 const controls=cells.map((row,k)=>{
  const laneOffset=mul(axes.u,(k/(N-1)-.5)*(uHi-uLo));
  const station=row.map((c,j)=>c.w>1e-8?mul(c.s,1/c.w):add(gc[j],mul(laneOffset,1-j/3)));
  const [a,b,c,d]=station,x=sub(mul(b,27),add(mul(a,8),d)),y=sub(mul(c,27),add(a,mul(d,8)));
  let p1=mul(sub(mul(x,2),y),1/18),p2=mul(sub(mul(y,2),x),1/18);const chord=len(sub(d,a)),cap=Math.max(chord*.55,.02);
  const bound=(p,base)=>{const q=sub(p,base),m=len(q);return add(base,mul(q,Math.min(1,cap/Math.max(m,1e-9))))};
  return [a,bound(p1,mix(a,d,1/3)),bound(p2,mix(a,d,2/3)),d];
 });
 if(controls.flat(2).some(x=>!Number.isFinite(x)))return null;
 const refs=[],seen=new Set(),triCount=faces.length/3,step=Math.max(1,Math.ceil(triCount/80));
 for(let f=0;f<triCount;f+=step){const ids=faces.slice(f*3,f*3+3);for(const [i,j] of [[ids[0],ids[1]],[ids[1],ids[2]],[ids[2],ids[0]]]){const key=Math.min(i,j)+':'+Math.max(i,j);if(seen.has(key))continue;seen.add(key);const dist=len(sub(vertices[i],vertices[j]));if(dist>=.002)refs.push({p:vertices[i].slice(),q:vertices[j].slice(),w:[weights[i],weights[j]],l:[uv[i],uv[j]],d:dist})}}
 const tris=[],triStep=Math.max(1,Math.ceil(triCount/56));
 for(let f=Math.floor(triStep/2);f<triCount;f+=triStep){const ids=faces.slice(f*3,f*3+3),p=ids.map(i=>vertices[i].slice()),area=len(cross(sub(p[1],p[0]),sub(p[2],p[0])));if(area>=.000004)tris.push({p,w:ids.map(i=>weights[i]),l:ids.map(i=>uv[i]),area})}
 const bounds=[],stride=Math.max(1,Math.ceil(vertices.length/220));for(let i=0;i<vertices.length;i+=stride)bounds.push({p:vertices[i].slice(),w:weights[i],l:uv[i]});
 return {...field,controls,uv,refs,tris,bounds};
}
function pose(field,a,b,wrap){
 const rows=field.controls.map(c=>c.map((p,j)=>{if(j===0)return apply(a,p);if(j===3)return apply(b,p);const t=j/3,p1=mix(apply(a,p),apply(b,p),smooth(t));return wrap?wrap(p,p1,t,a,b):p1}));
 const out=rows.map(r=>r.map(p=>p.slice()));
 for(let j=1;j<=2;j++)for(let k=0;k<N;k++){const ds=[];for(let q=Math.max(0,k-1);q<=Math.min(N-1,k+1);q++)ds.push(sub(rows[q][j],field.controls[q][j]));let avg=[0,0,0];for(const d of ds)avg=add(avg,d);avg=mul(avg,1/ds.length);out[k][j]=add(rows[k][j],mul(sub(avg,sub(rows[k][j],field.controls[k][j])),.10))}
 return out;
}
function frame(controls,l,t,reference=null){
 const c=laneCurve(controls,l),center=curve(c,t),T=tangent(c,t),dl=.25,l0=clamp(l-dl,0,N-1),l1=clamp(l+dl,0,N-1),U0=curve(laneCurve(controls,l0),t),U1=curve(laneCurve(controls,l1),t),U=sub(U1,U0);
 if(len(T)<1e-7||len(U)<1e-7)return null;let tn=norm(T);if(reference?.t&&dot(tn,reference.t)<0)tn=mul(tn,-1);let un=sub(U,mul(tn,dot(U,tn)));if(len(un)<1e-7)return null;un=norm(un);if(reference?.u&&dot(un,reference.u)<0)un=mul(un,-1);let nn=norm(cross(tn,un));if(len(nn)<.5)return null;if(reference?.n&&dot(nn,reference.n)<0){un=mul(un,-1);nn=mul(nn,-1)}
 return {center,T,U,t:tn,u:un,n:nn,tLen:len(T),uLen:len(U)/Math.max(l1-l0,1e-6)};
}
function point(p,t,l,field,posed,a,b,base,alpha=1){
 if(t<=0||t>=1||alpha<=0)return base;const R=frame(field.controls,l,t,{n:field.n});if(!R)return base;const q=qmix(a.q,b.q,smooth(t)),S=frame(posed,l,t,{t:qrot(q,R.t),u:qrot(q,R.u),n:qrot(q,R.n)});if(!S)return base;
 const off=sub(p,R.center),xt=dot(off,R.t),xu=dot(off,R.u),xn=dot(off,R.n),st=clamp(S.tLen/Math.max(R.tLen,1e-8),.68,1.45),su=clamp(S.uLen/Math.max(R.uLen,1e-8),.72,1.30),sn=clamp(1/(st*su),.68,1.35);
 const candidate=add(S.center,add(add(mul(S.t,xt*st),mul(S.u,xu*su)),mul(S.n,xn*sn))),gate=smooth(t/.14)*smooth((1-t)/.16)*clamp(alpha,0,1);
 return mix(base,candidate,gate);
}
function lane(p,t,field){return laneValue(p,t,field)}
function assess(field,posed,a,b,baseline,constraint=null){
 const rows=field.refs.map(e=>{const p=baseline(e.p,e.w[0]),q=baseline(e.q,e.w[1]);return {p,q,P:point(e.p,e.w[0],e.l[0],field,posed,a,b,p,1),Q:point(e.q,e.w[1],e.l[1],field,posed,a,b,q,1),d:e.d}});
 const faces=(field.tris||[]).map(e=>{const p=e.p.map((x,i)=>baseline(x,e.w[i]));return {p,q:e.p.map((x,i)=>point(x,e.w[i],e.l[i],field,posed,a,b,p[i],1)),rest:e.p,area:e.area}});
 const contacts=constraint?field.bounds.filter((_,i)=>i%3===0).map(e=>{const p=baseline(e.p,e.w);return {p,q:point(e.p,e.w,e.l,field,posed,a,b,p,1),rest:e.p}}):[];
 const penalty=r=>{const x=Math.abs(Math.log(Math.max(r,1e-6)));return x*x+4*Math.max(0,x-Math.log(1.9))**2};
 const score=alpha=>{let sum=0,max=0,areaSum=0,areaMax=0,contact=0,nc=0,minRatio=Infinity,maxRatio=0;
  for(const e of rows){const r=len(sub(mix(e.p,e.P,alpha),mix(e.q,e.Q,alpha)))/e.d;if(!Number.isFinite(r))return {loss:Infinity,areaLoss:Infinity,contactLoss:Infinity,minRatio:0,maxRatio:Infinity};const p=penalty(r);sum+=p;max=Math.max(max,p);minRatio=Math.min(minRatio,r);maxRatio=Math.max(maxRatio,r)}
  for(const e of faces){const p=e.p.map((x,i)=>mix(x,e.q[i],alpha)),r=len(cross(sub(p[1],p[0]),sub(p[2],p[0])))/e.area;if(!Number.isFinite(r))return {loss:Infinity,areaLoss:Infinity,contactLoss:Infinity,minRatio:0,maxRatio:Infinity};const q=penalty(r);areaSum+=q;areaMax=Math.max(areaMax,q)}
  for(const e of contacts){const v=constraint(e.rest,mix(e.p,e.q,alpha));if(Number.isFinite(v)){contact+=v*v;nc++}}
  return {loss:rows.length?sum/rows.length+.08*max:0,areaLoss:faces.length?areaSum/faces.length+.08*areaMax:null,contactLoss:nc?contact/nc:null,minRatio,maxRatio};
 };
 const base=score(0),safe=c=>Number.isFinite(c.loss)&&c.loss<=base.loss+.035&&(base.areaLoss===null||c.areaLoss!==null&&c.areaLoss<=base.areaLoss+.04)&&(base.contactLoss===null||c.contactLoss!==null&&c.contactLoss<=base.contactLoss+1e-8);
 let alpha=0,best=base;for(const x of [1,.85,.7,.55,.4,.25]){const c=score(x);if(safe(c)){alpha=x;best=c;break}}
 return {alpha,baseline:base.loss,score:best.loss,minRatio:best.minRatio,maxRatio:best.maxRatio,samples:rows.length,areaSamples:faces.length,areaBaseline:base.areaLoss,areaScore:best.areaLoss,contactBaseline:base.contactLoss,contactScore:best.contactLoss};
}
root.PectoralFanEngineV070=Object.freeze({version:'0.70.0',N,build,pose,point,lane,assess,curve,tangent,laneCurve,frame});
})(globalThis);
