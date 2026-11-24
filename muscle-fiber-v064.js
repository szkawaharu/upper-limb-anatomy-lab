/* FiberEngineV064 — deterministic geometric surface deformation.
 * No forces, activation, collision solver, measured fascicles or clinical validation.
 * Six continuous regional guide lanes, cubic centerlines, bounded local stretch.
 * The atlas rest surface is never replaced. CPU reference also drives audit.
 */
(function(root){
'use strict';
const N=6,clamp=(x,a,b)=>Math.max(a,Math.min(b,x)),smooth=x=>{x=clamp(x,0,1);return x*x*(3-2*x)};
const add=(a,b)=>a.map((x,i)=>x+b[i]),sub=(a,b)=>a.map((x,i)=>x-b[i]),mul=(a,s)=>a.map(x=>x*s),dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0),cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],len=a=>Math.hypot(...a),norm=a=>mul(a,1/Math.max(len(a),1e-12)),mix=(a,b,t)=>add(mul(a,1-t),mul(b,t));
const qrot=(q,v)=>{const u=q.slice(0,3),t=mul(cross(u,v),2);return add(v,add(mul(t,q[3]),cross(u,t)))};
const qmix=(a,b,t)=>norm(mix(a,mul(b,dot(a,b)<0?-1:1),t));
const apply=(a,p)=>add(a.p,qrot(a.q,sub(p,a.rest)));
function align(u,v,p){
 const d=clamp(dot(u,v),-1,1);
 if(d<-.995)return p.slice(); // reject ill-defined 180-degree correction; audit will decide
 const k=cross(u,v),a=cross(k,p);
 return add(p,mul(add(a,mul(cross(k,a),1/(1+d))),smooth((d+.995)/.045)));
}
function curve(c,t){const u=1-t;return add(add(mul(c[0],u*u*u),mul(c[1],3*u*u*t)),add(mul(c[2],3*u*t*t),mul(c[3],t*t*t)))}
function tangent(c,t){const u=1-t;return add(add(mul(sub(c[1],c[0]),3*u*u),mul(sub(c[2],c[1]),6*u*t)),mul(sub(c[3],c[2]),3*t*t))}
function lane(p,field){return clamp((p[field.axis]-field.lo)/Math.max(field.hi-field.lo,1e-8),0,1)*(N-1)}
function laneCurve(controls,l){const i=Math.min(N-2,Math.floor(clamp(l,0,N-1))),f=clamp(l-i,0,1);return controls[i].map((p,j)=>mix(p,controls[i+1][j],f))}
function build(vertices,faces,weights,axis){
 if(![0,1,2].includes(axis)||vertices.length!==weights.length||vertices.length<8||!faces.length)return null;
 let lo=Infinity,hi=-Infinity;
 for(const p of vertices){lo=Math.min(lo,p[axis]);hi=Math.max(hi,p[axis])}
 if(hi-lo<1e-7)return null;
 const field={axis,lo,hi},cells=Array.from({length:N},()=>Array.from({length:4},()=>({s:[0,0,0],n:0}))),global=Array.from({length:4},()=>({s:[0,0,0],n:0})),uv=new Float32Array(vertices.length);
 const addCell=(c,p,w)=>{if(w<=0)return;c.s=add(c.s,mul(p,w));c.n+=w};
 for(let i=0;i<vertices.length;i++){
  const p=vertices[i],l=lane(p,field),u=clamp(weights[i],0,1)*3;uv[i]=l;
  for(let j=Math.max(0,Math.floor(u)-1);j<=Math.min(3,Math.ceil(u)+1);j++){
   const w=Math.max(0,1-Math.abs(j-u)/1.15)**2;
   addCell(global[j],p,w);
   for(let k=Math.max(0,Math.floor(l)-1);k<=Math.min(N-1,Math.ceil(l)+1);k++)addCell(cells[k][j],p,w*Math.max(0,1-Math.abs(k-l)/1.4)**2);
  }
 }
 if(global.some(c=>c.n<1e-8))return null;
 const gc=global.map(c=>mul(c.s,1/c.n));
 const controls=cells.map(row=>{
  let offset=[0,0,0],n=0;for(let j=0;j<4;j++)if(row[j].n>1e-8){offset=add(offset,sub(mul(row[j].s,1/row[j].n),gc[j]));n++}
  offset=mul(offset,1/Math.max(1,n));
  const stations=row.map((c,j)=>c.n>1e-8?mul(c.s,1/c.n):add(gc[j],offset));
  // Interpolate the regional centroids at t=1/3 and 2/3, with a moderate
  // tangent limiter to avoid severe overshoot from sparse atlas sampling.
  const [a,b,c,d]=stations, x=sub(mul(b,27),add(mul(a,8),d)),y=sub(mul(c,27),add(a,mul(d,8)));
  const p1=mul(sub(mul(x,2),y),1/18),p2=mul(sub(mul(y,2),x),1/18);
  const chord=len(sub(d,a)),cap=Math.max(chord*.7,.02);
  function bound(p,base){const v=sub(p,base),n=len(v);return add(base,mul(v,Math.min(1,cap/Math.max(n,1e-8))))}
  return [a,bound(p1,mix(a,d,1/3)),bound(p2,mix(a,d,2/3)),d];
 });
 if(controls.some(c=>len(sub(c[3],c[0]))<.01||c.flat().some(x=>!Number.isFinite(x))))return null;
 const refs=[],seen=new Set(),step=Math.max(1,Math.ceil(faces.length/3/60));
 for(let i=0;i<faces.length/3;i+=step)for(const [j,k] of [[faces[i*3],faces[i*3+1]],[faces[i*3+1],faces[i*3+2]],[faces[i*3+2],faces[i*3]]]){
  const key=Math.min(j,k)+':'+Math.max(j,k);if(seen.has(key))continue;seen.add(key);
  const p=vertices[j],q=vertices[k],d=len(sub(p,q));if(d<.002)continue;
  refs.push({p:p.slice(),q:q.slice(),w:[weights[j],weights[k]],l:[uv[j],uv[k]],d});
 }
 const bounds=[],stride=Math.max(1,Math.ceil(vertices.length/180));
 for(let i=0;i<vertices.length;i+=stride)bounds.push({p:vertices[i].slice(),w:weights[i],l:uv[i]});
 // Explicitly retain coordinate extrema as well as regular samples.
 for(let axis=0;axis<3;axis++)for(const sign of [-1,1]){
  let j=0;for(let i=1;i<vertices.length;i++)if(sign*vertices[i][axis]>sign*vertices[j][axis])j=i;
  bounds.push({p:vertices[j].slice(),w:weights[j],l:uv[j]});
 }
 return {...field,controls,uv,refs,bounds};
}
function pose(field,a,b,wrap){
 return field.controls.map(c=>{
  const A=apply(a,c[0]),B=apply(b,c[3]),r=sub(c[3],c[0]),s=sub(B,A);
  return c.map((p,j)=>{
   if(j===0)return A;if(j===3)return B;
   const t=j/3,q=qmix(a.q,b.q,t),R=norm(qrot(q,r)),S=norm(s);
   let p1=len(s)>.01?add(mix(A,B,t),align(R,S,qrot(q,sub(p,mix(c[0],c[3],t))))):mix(apply(a,p),apply(b,p),t);
   if(wrap)p1=wrap(p,p1,t,a,b);
   return p1;
  });
 });
}
function point(p,t,l,field,posed,a,b,base,alpha=1){
 if(t<=0||t>=1||alpha<=0)return base;
 const r=laneCurve(field.controls,l),s=laneCurve(posed,l),u=tangent(r,t),v=tangent(s,t);
 const lr=len(u),ls=len(v);if(lr<1e-5||ls<1e-5)return base;
 const ru=mul(u,1/lr),sv=mul(v,1/ls),q=qmix(a.q,b.q,t),ref=qrot(q,ru);
 const offset=sub(p,curve(r,t)),long=dot(offset,ru),ax=clamp(ls/lr,.5,2),trans=1/Math.sqrt(ax);
 const scaled=add(mul(sub(offset,mul(ru,long)),trans),mul(ru,long*ax));
 const candidate=add(curve(s,t),align(ref,sv,qrot(q,scaled)));
 const gate=smooth(t/.18)*smooth((1-t)/.18)*smooth((dot(ref,sv)+.995)/.045);
 return mix(base,candidate,gate*alpha);
}
function assess(field,posed,a,b,baseline){
 const rows=field.refs.map(e=>{
  const p=baseline(e.p,e.w[0]),q=baseline(e.q,e.w[1]);
  return {p,q,P:point(e.p,e.w[0],e.l[0],field,posed,a,b,p),Q:point(e.q,e.w[1],e.l[1],field,posed,a,b,q),d:e.d};
 });
 const score=alpha=>{
  let sum=0,max=0,minRatio=Infinity,maxRatio=0;
  for(const e of rows){const r=len(sub(mix(e.P,e.p,1-alpha),mix(e.Q,e.q,1-alpha)))/e.d;if(!Number.isFinite(r))return {loss:Infinity,minRatio:0,maxRatio:Infinity};const x=Math.abs(Math.log(Math.max(r,1e-6))),pen=x*x+4*Math.max(0,x-Math.log(2))**2;sum+=pen;max=Math.max(max,pen);minRatio=Math.min(minRatio,r);maxRatio=Math.max(maxRatio,r)}
  return {loss:rows.length?sum/rows.length+.1*max:Infinity,minRatio,maxRatio};
 };
 const base=score(0);let best=base,alpha=0;
 for(const x of [.25,.5,.75,1]){const c=score(x);if(c.loss<best.loss-1e-10){best=c;alpha=x}}
 return {alpha,baseline:base.loss,score:best.loss,minRatio:best.minRatio,maxRatio:best.maxRatio,samples:rows.length};
}
function tubeVolume(pts){let v=0;for(let i=1;i<pts.length;i++){const a=pts[i-1],b=pts[i];v+=Math.PI*len(sub(b.slice(0,3),a.slice(0,3)))*(a[3]**2+a[3]*b[3]+b[3]**2)/3}return v}
function tubeCorrect(pts,target){
 if(!Number.isFinite(target)||target<=0||pts.length<3)return pts;
 const scale=k=>pts.map((p,i)=>[...p.slice(0,3),p[3]*(i===0||i===pts.length-1?1:k)]);
 let lo=.8,hi=1.25;for(let j=0;j<12;j++){const m=(lo+hi)/2;if(tubeVolume(scale(m))<target)lo=m;else hi=m}
 return scale((lo+hi)/2);
}
root.FiberEngineV064=Object.freeze({version:'0.64.0',N,build,lane,pose,point,assess,tubeVolume,tubeCorrect,curve,tangent});
})(globalThis);
