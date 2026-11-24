/* SerratusEngineV065 — regional serratus anterior slip deformation.
 * Ribs 1–9 are represented as nine continuous guide slips.  This is an
 * anatomical-region model, not segmented fascicles, force/activation, FE,
 * or measured subject-specific muscle architecture.
 */
(function(root){
'use strict';
const N=9,clamp=(x,a,b)=>Math.max(a,Math.min(b,x)),smooth=x=>{x=clamp(x,0,1);return x*x*(3-2*x)};
const add=(a,b)=>a.map((x,i)=>x+b[i]),sub=(a,b)=>a.map((x,i)=>x-b[i]),mul=(a,s)=>a.map(x=>x*s),dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0),cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],len=a=>Math.hypot(...a),norm=a=>mul(a,1/Math.max(len(a),1e-12)),mix=(a,b,t)=>add(mul(a,1-t),mul(b,t));
const qrot=(q,v)=>{const u=q.slice(0,3),t=mul(cross(u,v),2);return add(v,add(mul(t,q[3]),cross(u,t)))};
const qmix=(a,b,t)=>norm(mix(a,mul(b,dot(a,b)<0?-1:1),t));
const apply=(a,p)=>add(a.p,qrot(a.q,sub(p,a.rest)));
function align(u,v,p){const d=clamp(dot(u,v),-1,1);if(d<-.995)return p.slice();const k=cross(u,v),a=cross(k,p);return add(p,mul(add(a,mul(cross(k,a),1/(1+d))),smooth((d+.995)/.045)))}
function curve(c,t){const u=1-t;return add(add(mul(c[0],u*u*u),mul(c[1],3*u*u*t)),add(mul(c[2],3*u*t*t),mul(c[3],t*t*t)))}
function tangent(c,t){const u=1-t;return add(add(mul(sub(c[1],c[0]),3*u*u),mul(sub(c[2],c[1]),6*u*t)),mul(sub(c[3],c[2]),3*t*t))}
function polyline(path,t){
 if(path.length===1)return path[0].slice();const seg=[],cum=[0];let total=0;
 for(let i=1;i<path.length;i++){const d=len(sub(path[i],path[i-1]));seg.push(d);total+=d;cum.push(total)}
 const s=clamp(t,0,1)*Math.max(total,1e-9);let i=0;while(i<seg.length-1&&cum[i+1]<s)i++;
 return mix(path[i],path[i+1],seg[i]?clamp((s-cum[i])/seg[i],0,1):0);
}
function sampled(list,limit=180){const step=Math.max(1,Math.ceil(list.length/limit));return list.filter((_,i)=>i%step===0)}
function originAnchors(vertices,weights,ribs){
 const low=sampled(vertices.filter((_,i)=>weights[i]<=.30),420);if(low.length<8)return null;
 const out=[];
 for(const rib of ribs){
  const pts=sampled(rib.points||[],180);if(!pts.length)return null;
  let best=null,d=Infinity;
  for(const q of pts)for(const p of low){const x=len(sub(p,q));if(x<d){d=x;best=q}}
  if(!best||!Number.isFinite(d))return null;
  out.push({level:rib.level,p:best.slice(),distance:d});
 }
 return out.sort((a,b)=>a.level-b.level);
}
function slipByY(y,origins){
 // Rib anchors descend with increasing level. Interpolate continuously so
 // adjacent slips share the surface instead of creating visible seams.
 const ys=origins.map(x=>x.p[1]);if(y>=ys[0])return 0;if(y<=ys[N-1])return N-1;
 for(let i=0;i<N-1;i++)if(y<=ys[i]&&y>=ys[i+1])return i+clamp((ys[i]-y)/Math.max(ys[i]-ys[i+1],1e-9),0,1);
 let best=0;for(let i=1;i<N;i++)if(Math.abs(y-ys[i])<Math.abs(y-ys[best]))best=i;return best;
}
function insertionT(level){
 // Literature supports upper (1–2), middle (3–6) and lower (6–9)
 // functional regions; values below interpolate those regions, not measured
 // individual fascicle insertion coordinates.
 return [0,.07,.22,.37,.51,.65,.78,.90,1][level-1];
}
function build(vertices,faces,weights,ribs,medialPath){
 if(!vertices?.length||vertices.length!==weights?.length||!faces?.length||!medialPath?.length||ribs?.length!==N)return null;
 const origins=originAnchors(vertices,weights,ribs);if(!origins)return null;
 const slip=new Float32Array(vertices.length);for(let i=0;i<vertices.length;i++)slip[i]=slipByY(vertices[i][1],origins);
 const controls=[];
 for(let k=0;k<N;k++){
  const a=origins[k].p.slice(),d=polyline(medialPath,insertionT(k+1));
  const c=[a,null,null,d],sum=[[0,0,0],[0,0,0]],den=[0,0];
  for(let i=0;i<vertices.length;i++){
   const lane=Math.max(0,1-Math.abs(slip[i]-k)/1.1),w=weights[i];if(lane<=0)continue;
   for(let j=0;j<2;j++){const t=j?2/3:1/3,q=Math.max(0,1-Math.abs(w-t)/.24)**2*lane;if(q){sum[j]=add(sum[j],mul(vertices[i],q));den[j]+=q}}
  }
  c[1]=den[0]?mul(sum[0],1/den[0]):mix(a,d,1/3);c[2]=den[1]?mul(sum[1],1/den[1]):mix(a,d,2/3);
  const chord=len(sub(d,a)),cap=Math.max(.04,chord*.65);
  for(let j=1;j<=2;j++){const base=mix(a,d,j/3),v=sub(c[j],base),n=len(v);if(n>cap)c[j]=add(base,mul(v,cap/n))}
  controls.push(c);
 }
 const refs=[],seen=new Set(),step=Math.max(1,Math.ceil(faces.length/3/220));
 for(let t=0;t<faces.length/3&&refs.length<220;t+=step){
  const tri=[faces[3*t],faces[3*t+1],faces[3*t+2]];
  for(const [i,j] of [[tri[0],tri[1]],[tri[1],tri[2]],[tri[2],tri[0]]]){
   const key=Math.min(i,j)+':'+Math.max(i,j);if(seen.has(key))continue;seen.add(key);
   const d=len(sub(vertices[i],vertices[j]));if(d<.001)continue;
   refs.push({p:vertices[i].slice(),q:vertices[j].slice(),w:[weights[i],weights[j]],s:[slip[i],slip[j]],d});
   if(refs.length>=220)break;
  }
 }
 const bounds=[],stride=Math.max(1,Math.ceil(vertices.length/220));for(let i=0;i<vertices.length;i+=stride)bounds.push({p:vertices[i].slice(),w:weights[i],s:slip[i]});
 return {controls,slip,origins,refs,bounds,medialPath:medialPath.map(x=>x.slice())};
}
function laneCurve(controls,l){const i=Math.min(N-2,Math.floor(clamp(l,0,N-1))),f=clamp(l-i,0,1);return controls[i].map((p,j)=>mix(p,controls[i+1][j],f))}
function pose(field,trunk,scapula,wrap){
 return field.controls.map(c=>c.map((p,j)=>{
  const t=j/3;if(j===0)return apply(trunk,p);if(j===3)return apply(scapula,p);
  let q=mix(apply(trunk,p),apply(scapula,p),smooth(t));
  if(wrap)q=wrap(p,q,t);return q;
 }));
}
function point(p,t,s,field,posed,trunk,scapula,base,alpha=1){
 if(t<=0||t>=1||alpha<=0)return base;
 const r=laneCurve(field.controls,s),n=laneCurve(posed,s),u=tangent(r,t),v=tangent(n,t),lr=len(u),lv=len(v);if(lr<1e-6||lv<1e-6)return base;
 const ru=mul(u,1/lr),sv=mul(v,1/lv),q=qmix(trunk.q,scapula.q,smooth(t)),ref=norm(qrot(q,ru));
 const off=sub(p,curve(r,t)),long=dot(off,ru),ax=clamp(lv/lr,.55,1.8),trans=1/Math.sqrt(ax);
 const scaled=add(mul(sub(off,mul(ru,long)),trans),mul(ru,long*ax)),candidate=add(curve(n,t),align(ref,sv,qrot(q,scaled)));
 const gate=smooth(t/.16)*smooth((1-t)/.16)*smooth((dot(ref,sv)+.995)/.045);
 return mix(base,candidate,gate*alpha);
}
function assess(field,posed,trunk,scapula,baseline,constraint=null){
 const rows=field.refs.map(e=>{const p=baseline(e.p,e.w[0]),q=baseline(e.q,e.w[1]);return {...e,P:point(e.p,e.w[0],e.s[0],field,posed,trunk,scapula,p),Q:point(e.q,e.w[1],e.s[1],field,posed,trunk,scapula,q),bp:p,bq:q}});
 const score=a=>{let loss=0,max=0,contact=0,n=0,min=Infinity,hi=0;for(const e of rows){const P=mix(e.bp,e.P,a),Q=mix(e.bq,e.Q,a),r=len(sub(P,Q))/e.d,x=Math.abs(Math.log(Math.max(r,1e-6))),z=x*x+4*Math.max(0,x-Math.log(1.8))**2;loss+=z;max=Math.max(max,z);min=Math.min(min,r);hi=Math.max(hi,r);if(constraint){for(const [rest,p] of [[e.p,P],[e.q,Q]]){const c=constraint(rest,p);if(Number.isFinite(c)){contact+=c*c;n++}}}}return {loss:rows.length?loss/rows.length+.15*max:Infinity,contact:n?contact/n:null,minRatio:min,maxRatio:hi}};
 const base=score(0);let best=base,alpha=0;
 for(const a of [.2,.4,.6,.8,1]){const c=score(a),safe=c.loss<=base.loss+1e-10&&(base.contact===null||c.contact!==null&&c.contact<=base.contact+1e-8);if(safe&&c.loss<best.loss-1e-10){best=c;alpha=a}}
 return {alpha,baseline:base.loss,score:best.loss,contactBaseline:base.contact,contactScore:best.contact,minRatio:best.minRatio,maxRatio:best.maxRatio,samples:rows.length};
}
root.SerratusEngineV065=Object.freeze({version:'0.65.0',N,build,pose,point,assess,laneCurve,insertionT});
})(globalThis);
