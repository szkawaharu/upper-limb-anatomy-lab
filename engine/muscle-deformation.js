(function(root){'use strict';const M=root.AnatomyMath,P=root.FiberPath,{add,sub,mul,dot,len,mix,qrot,apply,smooth}=M;
function bind(definition,weights,regionalField=null){
 const {vertices,faces}=definition.restMesh;
 if(vertices.length!==weights.length||vertices.some(p=>!p.every(Number.isFinite))||weights.some(w=>!Number.isFinite(w)||w<0||w>1))throw Error('Invalid muscle data '+definition.id);
 const sheet=['fan','sheet','pennate'].includes(definition.morphology),fan=sheet&&root.PectoralFanEngineV070?.build(vertices,faces,weights);
 const f=regionalField||fan||root.FiberEngineV064?.build(vertices,faces,weights,1);if(!f)throw Error('No regional field '+definition.id);
 const lane=f.slip||f.uv;if(!lane)throw Error('Missing material coordinates');
 const ref=fan?.u||f.u||[0,0,1];
 // Attachment blend weight is not a material arc coordinate. Find the
 // closest rest-centerline location once, avoiding a spurious longitudinal
 // residual that would be stretched into the sheet normal during elevation.
 const material=new Float32Array(weights.length);
 for(let i=0;i<weights.length;i++){const c=P.laneCurve(f.controls,lane[i]),p=vertices[i],lo=Math.max(0,weights[i]-.25),hi=Math.min(1,weights[i]+.25);let best=weights[i],loss=Infinity;for(let k=0;k<=24;k++){const t=lo+(hi-lo)*k/24,d=len(sub(p,P.curve(c,t)));if(d<loss){loss=d;best=t}}for(let k=0;k<5;k++){const q=P.curve(c,best),v=P.tangent(c,best),step=dot(sub(q,p),v)/Math.max(dot(v,v),1e-12);best=Math.max(lo,Math.min(hi,best-step))}material[i]=best}

 const bounds=new Set();for(let i=0;i<vertices.length;i+=Math.max(1,Math.floor(vertices.length/256)))bounds.add(i);for(let k=0;k<3;k++)for(const sign of [-1,1]){let j=0;for(let i=1;i<vertices.length;i++)if(sign*vertices[i][k]>sign*vertices[j][k])j=i;bounds.add(j)}
 return {bounds:[...bounds],definition,weights:Float32Array.from(weights),field:f,lane,material,sheet,reference:ref,rest:Float32Array.from(vertices.flat()),faces:Uint32Array.from(faces),cache:new Map()};
}
function pose(binding,a,b,contact=null,originAt=null){
 const f=binding.field,posed=f.controls.map(c=>c.map((p,j)=>mix(apply(originAt?originAt(p):a,p),apply(b,p),smooth(j/3))));
 return {binding,a,b,posed,contact,originAt,restFlat:Float32Array.from(f.controls.flat(2)),nowFlat:Float32Array.from(posed.flat(2)),activation:null,force:null};
}
function evaluate(state,p,t,l,normal=null,materialT=t){
 const {binding:B,b,posed}=state,a=state.originAt?state.originAt(p):state.a,base=mix(apply(a,p),apply(b,p),t),baseN=normal?mix(qrot(a.q,normal),qrot(b.q,normal),t):null;
 if(t<=0||t>=1)return {point:base,normal:baseN,stretch:1,determinant:1,bounded:false};
 const R=P.frame(B.field.controls,l,materialT,B.reference),raw=P.frame(posed,l,materialT,qrot(a.q,B.reference));
 if(!R||!raw)return {point:base,normal:baseN,invalidFrame:true,stretch:null,determinant:null};
 const S=P.transported(R,raw,a,b,materialT),tendon=root.TendonDeformation.fraction(t,B.definition.tissues),sc=root.VolumeConstraint.scales(S.tl/R.tl,raw.ul/Math.max(R.ul,1e-7),B.sheet,tendon),off=sub(p,R.c),coords=[dot(off,R.t),dot(off,R.u),dot(off,R.n)],candidate=add(S.c,add(add(mul(S.t,coords[0]*sc.t),mul(S.u,coords[1]*sc.u)),mul(S.n,coords[2]*sc.n))),gate=smooth(t/.12)*smooth((1-t)/.12);
 let point=mix(base,candidate,gate);point=root.ContactConstraint.project(p,point,t,state.contact);
 const n=normal?mix(baseN,add(add(mul(S.t,dot(normal,R.t)/sc.t),mul(S.u,dot(normal,R.u)/sc.u)),mul(S.n,dot(normal,R.n)/sc.n)),gate):null;
 if(!point.every(Number.isFinite))throw Error('Non-finite deformation '+B.definition.id);
 return {point,normal:n,stretch:sc.t,determinant:sc.determinant,bounded:sc.bounded};
}
function point(state,p,t,l,materialT=t){return evaluate(state,p,t,l,null,materialT).point}
function deform(state,output=new Float32Array(state.binding.rest.length)){
 const B=state.binding;if(output.length!==B.rest.length)throw Error('Output length');
 for(let i=0;i<B.weights.length;i++)output.set(point(state,B.rest.subarray(3*i,3*i+3),B.weights[i],B.lane[i],B.material[i]),3*i);return output;
}
// A registry exposes declared adapters honestly; procedural route geometry is
// preserved until a measured rest mesh/attachment field can be supplied.
class Engine {
 constructor(){this.bindings=new Map();this.adapters=new Map();this.fallbacks=new Map()}
 register(binding){this.bindings.set(binding.definition.id,binding);return binding}
 registerAdapter(id,kind,metadata={}){this.adapters.set(id,{id,kind,...metadata,validation:'legacy-route-adapter',activeForce:false})}
 fallback(id,reason){const key=id+':'+reason;this.fallbacks.set(key,(this.fallbacks.get(key)||0)+1)}
 report(){return {version:'0.73.0',atlas:[...this.bindings.keys()],adapters:[...this.adapters.values()],fallbacks:[...this.fallbacks].map(([reason,count])=>({reason,count})),activation:null,force:null}}
}
root.MuscleDeformation=Object.freeze({version:'0.73.0',bind,pose,evaluate,point,deform,Engine});
})(globalThis);
