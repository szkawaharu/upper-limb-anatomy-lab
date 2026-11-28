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
root.GeodesicAttachmentField=Object.freeze({field,rebuild});
})(globalThis);
