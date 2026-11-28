/* Experimental mesh-geodesic coordinate candidate. No spatial nearest-neighbor
 * shortcuts, no change to the source positions/faces or fixed attachment masks.
 * Geodesic distances are edge-path upper bounds, not fascicle lengths. */
(function(root){
'use strict';
function field(vertices,faces,origin,insertion){
 const n=vertices.length,adj=Array.from({length:n},()=>new Map());
 if(!origin.length||!insertion.length)throw Error('Both attachment regions required');
 const pin=new Int8Array(n).fill(-1);
 for(const [k,ids] of [[0,origin],[1,insertion]])for(const i of ids){
  if(!Number.isInteger(i)||i<0||i>=n||pin[i]>=0&&pin[i]!==k)throw Error('Overlapping or invalid attachment mask');
  pin[i]=k;
 }
 for(let i=0;i<faces.length;i+=3)for(let j=0;j<3;j++){
  const a=faces[i+j],b=faces[i+(j+1)%3];
  if(a===b)continue;
  const d=Math.hypot(...vertices[a].map((x,k)=>x-vertices[b][k]));
  if(!(d>0)||!Number.isFinite(d))throw Error('Invalid source edge');
  adj[a].set(b,d);adj[b].set(a,d);
 }
 function distance(ids){
  const d=new Float64Array(n).fill(Infinity),heap=[];
  const push=(v,x)=>{let j=heap.length;heap.push([v,x]);while(j){const p=(j-1)>>1;if(heap[p][1]<=x)break;heap[j]=heap[p];j=p;}heap[j]=[v,x]};
  const pop=()=>{const a=heap[0],last=heap.pop();if(heap.length){let i=0;while(i*2+1<heap.length){let j=i*2+1;if(j+1<heap.length&&heap[j+1][1]<heap[j][1])j++;if(heap[j][1]>=last[1])break;heap[i]=heap[j];i=j;}heap[i]=last;}return a};
  for(const i of ids){d[i]=0;push(i,0)}
  while(heap.length){const [i,x]=pop();if(x!==d[i])continue;for(const [j,w] of adj[i])if(x+w<d[j]){d[j]=x+w;push(j,d[j])}}
  return d;
 }
 const a=distance(origin),b=distance(insertion),weights=new Float32Array(n);
 let minSpan=Infinity,maxGradient=0,fixedOppositeEdges=0;
 for(let i=0;i<n;i++){
  if(!Number.isFinite(a[i]+b[i])||a[i]+b[i]<=0)throw Error('Disconnected or conflicting attachment domain at '+i);
  minSpan=Math.min(minSpan,a[i]+b[i]);weights[i]=pin[i]>=0?pin[i]:a[i]/(a[i]+b[i]);
 }
 for(let i=0;i<n;i++)for(const [j,d] of adj[i])if(j>i){
  maxGradient=Math.max(maxGradient,Math.abs(weights[i]-weights[j])/d);
  if(pin[i]>=0&&pin[j]>=0&&pin[i]!==pin[j])fixedOppositeEdges++;
 }
 return {weights,minSpan,maxGradient,fixedOppositeEdges,physicalVolumeVerified:false};
}
function rebuild(B){
 const vertices=Array.from({length:B.weights.length},(_,i)=>Array.from(B.rest.subarray(3*i,3*i+3)));
 const original=B.definition.attachments;
 const g=field(vertices,B.faces,original[0].indices,original[1].indices);
 const def={...B.definition,restMesh:{vertices,faces:Array.from(B.faces)}};
 const next=root.MuscleDeformation.bind(def,g.weights);
 if(next.rest.length!==B.rest.length||next.faces.length!==B.faces.length||next.rest.some((x,i)=>x!==B.rest[i])||next.faces.some((x,i)=>x!==B.faces[i]))throw Error('Source changed');
 for(let i=0;i<B.weights.length;i++)if((B.weights[i]===0||B.weights[i]===1)&&next.weights[i]!==B.weights[i])throw Error('Existing exact pin changed');
 return {binding:next,diagnostics:{minSpanMm:g.minSpan/B.definition.source.worldUnitsPerMm,maxGradientPerMm:g.maxGradient*B.definition.source.worldUnitsPerMm,fixedOppositeEdges:g.fixedOppositeEdges,fieldMethod:next.fieldMethod,productAccepted:false}};
}
// Explicit anatomical-mask revision experiment, distinct from rebuild(), which
// preserves legacy pins. Coordinates and topology must match the audited mesh.
function rebindSurface(B,data){
 if(data.rest.length!==B.rest.length||data.faces.length!==B.faces.length||
    data.rest.some((x,i)=>x!==B.rest[i])||data.faces.some((x,i)=>x!==B.faces[i]))throw Error('Audited mesh mismatch');
 const vertices=Array.from({length:B.weights.length},(_,i)=>Array.from(B.rest.subarray(3*i,3*i+3)));
 const origin=data.origin,insertion=data.insertion;
 if(!origin.length||!insertion.length)throw Error('Missing audited attachment candidate');
 const pin=new Int8Array(vertices.length).fill(-1);
 for(const [k,ids] of [[0,origin],[1,insertion]])for(const i of ids){
  if(!Number.isInteger(i)||i<0||i>=vertices.length||pin[i]>=0)throw Error('Invalid audited mask');
  pin[i]=k;
 }
 let g,method='intrinsic-edge-distance',geodesicRejected=null;
 try{g=field(vertices,B.faces,origin,insertion)}
 catch(e){
  // Disconnected fragments cannot be interpolated along a connected path.
  // Explicit separate candidate using true nearest-surface distances, not
  // fictitious pins or discarded triangles.
  if(!String(e).includes('Disconnected or conflicting attachment domain'))throw e;
  geodesicRejected=String(e);method='source-triangle-distance-ratio';
  if(data.originDistances.length!==vertices.length||data.insertionDistances.length!==vertices.length)throw Error('Distance size mismatch');
  const weights=Float32Array.from(vertices,(_,i)=>{
   if(pin[i]>=0)return pin[i];
   const a=data.originDistances[i],b=data.insertionDistances[i];
   if(!(a>=0&&b>=0&&a+b>0&&Number.isFinite(a+b)))throw Error('Invalid source distances');
   return a/(a+b);
  });
  g={weights,minSpan:null,maxGradient:null,fixedOppositeEdges:null};
 }
 const attachments=[origin,insertion].map((ids,k)=>({...B.definition.attachments[k],indices:Uint32Array.from(ids),boneLocal:null,confidence:'source-near-surface-experimental'}));
 const def={...B.definition,attachments,attachmentFieldSource:{method,radiusMm:data.radiusMm,confidence:'near-surface candidate, not measured footprint'},restMesh:{vertices,faces:Array.from(B.faces)}};
 const next=root.MuscleDeformation.bind(def,g.weights);
 if(next.rest.some((x,i)=>x!==B.rest[i])||next.faces.some((x,i)=>x!==B.faces[i]))throw Error('Source changed');
 for(let i=0;i<pin.length;i++)if(pin[i]>=0&&next.weights[i]!==pin[i])throw Error('Revised pin changed');
 return {binding:next,diagnostics:{method,geodesicRejected,radiusMm:data.radiusMm,originCount:origin.length,insertionCount:insertion.length,originPerNearestC:data.originPerNearestC,oldPinCounts:B.definition.attachments.map(x=>x.indices.length),minSpanMm:g.minSpan===null?null:g.minSpan/B.definition.source.worldUnitsPerMm,maxGradientPerMm:g.maxGradient===null?null:g.maxGradient*B.definition.source.worldUnitsPerMm,fixedOppositeEdges:g.fixedOppositeEdges,fieldMethod:next.fieldMethod,productAccepted:false,neckMapping:'unchanged legacy mapping, not per-vertebra validated'}};
}
root.GeodesicAttachmentField=Object.freeze({field,rebuild,rebindSurface});
})(globalThis);
