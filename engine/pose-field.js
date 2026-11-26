/* Deterministic interpolation of the calibrated joint solution, never skin
 * vertex positions. FK still enforces SC/AC/GH connections exactly. */
(function(root){
'use strict';
function signature(value){
 const canon=x=>typeof x==='number'?Math.round(x*1e8)/1e8:Array.isArray(x)?x.map(canon):x&&typeof x==='object'?Object.fromEntries(Object.keys(x).sort().map(k=>[k,canon(x[k])])):x;
 const s=JSON.stringify(canon(value));let h=2166136261;
 for(let i=0;i<s.length;i++){h^=s.charCodeAt(i);h=Math.imul(h,16777619)}
 return 'fnv1a32-'+(h>>>0).toString(16).padStart(8,'0');
}
function sample(f,elev,plane){
 if(f.elevationKnots&&f.planeKnots)return nonuniform(f,elev,plane);
 const x=Math.max(0,Math.min(f.elevationMax,elev))/f.elevationStep,y=Math.max(0,Math.min(f.planeMax,plane))/f.planeStep,i=Math.min(f.elevationCount-2,Math.floor(x)),j=Math.min(f.planeCount-2,Math.floor(y)),u=x-i,v=y-j;
 const at=(e,p,c)=>f.values[(Math.max(0,Math.min(f.elevationCount-1,e))*f.planeCount+Math.max(0,Math.min(f.planeCount-1,p)))*6+c];
 const cubic=(a,b,c,d,t,first,last)=>{const m=first?0:(c-a)/2,n=last?0:(d-b)/2;return b+t*(m+t*(3*(c-b)-2*m-n+t*(2*(b-c)+m+n)))};
 const out=new Array(6);
 for(let c=0;c<6;c++){
  const r=[];for(let p=j-1;p<=j+2;p++)r.push(cubic(at(i-1,p,c),at(i,p,c),at(i+1,p,c),at(i+2,p,c),u,i===0,i===f.elevationCount-2));
  out[c]=cubic(...r,v,j===0,j===f.planeCount-2);
 }
 return out;
}
function nonuniform(f,elev,plane){
 const E=f.elevationKnots,P=f.planeKnots;
 const cell=(a,x)=>{let lo=0,hi=a.length-1;while(hi-lo>1){const m=(lo+hi)>>1;if(x<a[m])hi=m;else lo=m;}return Math.min(a.length-2,lo);};
 const x=Math.max(E[0],Math.min(E.at(-1),elev)),y=Math.max(P[0],Math.min(P.at(-1),plane)),i=cell(E,x),j=cell(P,y);
 const at=(e,p,c)=>f.values[(Math.max(0,Math.min(E.length-1,e))*P.length+Math.max(0,Math.min(P.length-1,p)))*6+c];
 const cubic=(a,b,c,d,knots,k,x)=>{const h=knots[k+1]-knots[k],t=(x-knots[k])/h,m=k===0?0:(c-a)*h/(knots[k+1]-knots[k-1]),n=k===knots.length-2?0:(d-b)*h/(knots[k+2]-knots[k]);return b+t*(m+t*(3*(c-b)-2*m-n+t*(2*(b-c)+m+n)));};
 return Array.from({length:6},(_,c)=>{const r=[];for(let p=j-1;p<=j+2;p++)r.push(cubic(at(i-1,p,c),at(i,p,c),at(i+1,p,c),at(i+2,p,c),E,i,x));return cubic(...r,P,j,y);});
}
root.ShoulderPoseField=Object.freeze({signature,sample,method:'joint-space-bicubic-nonuniform-v2'});
})(globalThis);
