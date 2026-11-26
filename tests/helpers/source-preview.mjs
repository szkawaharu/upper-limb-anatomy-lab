// Test-only hydration of captured source rest data. No refit or hidden muscles.
export const sourcePreview=String.raw`
const sourcePreviewRaw=MuscleDeformation;
const sourcePreviewPoint=(f,p,w,l,t)=>{
 const B=f.binding;if(!B.meshPreparation)return sourcePreviewRaw.point(f,p,w,l,t);const a=f.originAt?f.originAt(p):f.a,b=f.b;
 if(w<=0||w>=1)return AnatomyMath.apply(w>=1?b:a,p);
 const q=AnatomyMath.qblendLift(a.q,b.q,w),c=FiberPath.curve(FiberPath.laneCurve(B.field.controls,l),t),center=AnatomyMath.mix(AnatomyMath.apply(a,c),AnatomyMath.apply(b,c),w),raw=AnatomyMath.add(center,AnatomyMath.qrot(q,AnatomyMath.sub(p,c)));
 return ContactConstraint.project(p,raw,w,f.contact);
};
MuscleDeformation=Object.freeze({...sourcePreviewRaw,point:sourcePreviewPoint,deform(f){const B=f.binding;if(!B.meshPreparation)return sourcePreviewRaw.deform(f);const out=new Float32Array(B.rest.length);for(let i=0;i<B.weights.length;i++)out.set(sourcePreviewPoint(f,B.rest.subarray(i*3,i*3+3),B.weights[i],B.lane[i],B.material[i]),i*3);return out;}});
const sourcePreviewBinding=v73Binding;
v73Binding=(spec,parsed)=>{
 if(!parsed.prepared)return sourcePreviewBinding(spec,parsed);
 const x=parsed.prepared,a=v58Anchor(spec.a),b=v58Anchor(spec.b),definition=MuscleDefinitions.create(spec,parsed.data.vertices,parsed.data.faces,x.weights,[a,b]),B={...x,definition,rest:Float32Array.from(x.rest),faces:Uint32Array.from(x.faces),weights:Float32Array.from(x.weights),lane:Float32Array.from(x.lane),material:Float32Array.from(x.material),bounds:[],cache:new Map()};
 for(let i=0;i<B.weights.length;i+=Math.max(1,Math.floor(B.weights.length/256)))B.bounds.push(i);return v73System.register(B);
};
const sourcePreviewLoad=v58Load;
v58Load=async(fetcher=null)=>sourcePreviewLoad(fetcher||(async(spec)=>{
 const response=await fetch('/source-preview/'+spec.id+'.json.gz');if(!response.ok)throw Error('prepared source '+spec.id);
 const x=JSON.parse(await new Response(response.body.pipeThrough(new DecompressionStream('gzip'))).text()),vertices=Array.from({length:x.rest.length/3},(_,i)=>x.rest.slice(i*3,i*3+3)),data=mesh(vertices,x.faces),box=v58Box(vertices);
 return {prepared:x,data,box,weights:x.weights,commonWeights:x.weights,originalTriangles:x.meshPreparation.sourceTriangles,attachmentMode:'captured immutable source binding; no pose scoring'};
}));
$('v73Backend').value='cpu';
`;
