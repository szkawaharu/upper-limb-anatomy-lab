(function(root){'use strict';const {len,sub,clamp}=root.AnatomyMath;
class ProceduralMuscleAdapter {
 constructor(){this.rest=new Map();this.reports=new Map()}
 deform(id,path,neutral,key,ellipticity=1){const cacheKey=id+':'+key;let ref=this.rest.get(cacheKey);if(neutral&&!ref){ref=path.map(p=>p.slice());this.rest.set(cacheKey,ref)}if(!ref||ref.length!==path.length){this.reports.set(cacheKey,{id,status:'missing-neutral-reference',force:null});return path}
 const now=path.map(p=>p.slice());let bounded=0;
 for(let i=0;i<path.length;i++){const lo=Math.max(0,i-1),hi=Math.min(path.length-1,i+1),L0=len(sub(ref[hi],ref[lo])),L=len(sub(path[hi],path[lo]));if(i===0||i===path.length-1){now[i][3]=ref[i][3];continue}const raw=L/Math.max(L0,1e-8),s=clamp(1/Math.sqrt(clamp(raw,.35,2.8)),.6,1.6);if(raw<.35||raw>2.8)bounded++;now[i][3]=ref[i][3]*s}
 this.reports.set(cacheKey,{id,status:'passive-path-adapter',...root.TendonDeformation.metrics(ref,path),bounded,ellipticity});return now}
}
root.ProceduralMuscleAdapter=ProceduralMuscleAdapter;
})(globalThis);
