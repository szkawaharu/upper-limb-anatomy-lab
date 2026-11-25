/* Pose-independent graph regularization. Dirichlet attachment regions remain
 * exact. Neighbours come from mesh connectivity, never spatial proximity
 * across adjacent muscles or opposite sides of a thin sheet. */
(function(root){'use strict';
function adjacency(vertices,faces){const sets=Array.from({length:vertices.length},()=>new Set());for(let k=0;k<faces.length;k+=3)for(let j=0;j<3;j++){const a=faces[k+j],b=faces[k+(j+1)%3];sets[a].add(b);sets[b].add(a)}const offsets=new Uint32Array(vertices.length+1),indices=[],weights=[];for(let i=0;i<vertices.length;i++){offsets[i]=indices.length;let sum=0;const row=[...sets[i]];for(const j of row){const d=Math.hypot(...vertices[i].map((x,k)=>x-vertices[j][k])),w=1/Math.max(d,1e-6);indices.push(j);weights.push(w);sum+=w}for(let k=offsets[i];k<indices.length;k++)weights[k]/=sum||1}offsets[vertices.length]=indices.length;return {offsets,indices:Uint32Array.from(indices),weights:Float32Array.from(weights)}}
function regularize(values,attachmentWeights,graph,passes=80){const a=Float64Array.from(values),b=new Float64Array(a.length),{offsets,indices,weights}=graph;let current=a,next=b;for(let n=0;n<passes;n++){for(let i=0;i<a.length;i++){const w=attachmentWeights[i];if(w===0||w===1){next[i]=w;continue}let x=0;for(let k=offsets[i];k<offsets[i+1];k++)x+=current[indices[k]]*weights[k];next[i]=offsets[i+1]>offsets[i]?.99*x+.01*a[i]:a[i]}[current,next]=[next,current]}return Float32Array.from(current)}
root.MaterialCoordinates=Object.freeze({adjacency,regularize});
})(globalThis);
