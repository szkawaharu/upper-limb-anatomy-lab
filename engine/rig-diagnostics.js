(function(root){
'use strict';
const M=root.AnatomyMath,F=root.JointFrames,{add,sub,mul,norm,dot,cross,qrot,qmul,apply}=M;
function longBone(prox,dist,lateral){const y=norm(sub(prox,dist)),x=norm(cross(y,lateral)),z=norm(cross(x,y));return F.basis(x,y,z)}
function matrix(a){const x=qrot(a.q,[1,0,0]),y=qrot(a.q,[0,1,0]),z=qrot(a.q,[0,0,1]),t=sub(a.p,qrot(a.q,a.rest));return [...x,0,...y,0,...z,0,...t,1]}
function frames(model,bones){
 const L=model.landmarks,lateral=L.EL_LAT&&L.EL_MED?sub(L.EL_LAT,L.EL_MED):qrot(model.thoraxFrame,[0,0,1]),h=longBone(L.GH,L.EL,lateral),f=L.RH&&L.DU?longBone(L.RH,L.DU,lateral):h;
 const frame=(id,restQ,origin)=>{const a=bones[id],q=a?qmul(a.q,restQ):restQ;return {origin:a?apply(a,origin):origin,q,axes:[[1,0,0],[0,1,0],[0,0,1]].map(v=>qrot(q,v))}};
 return {thorax:frame('thorax',model.thoraxFrame,L.IJ||L.SC),clavicle:frame('clavicle',model.clavicleFrame,L.SC),scapula:frame('scapula',model.scapulaFrame,L.AA),humerus:frame('humerus',h,L.GH),...(bones.ulna?{ulna:frame('ulna',f,L.EL),radius:frame('radius',f,L.RH||L.EL)}:{})};
}
function connections(model,bones,radialHead=null,distalUlna=null){
 const L=model.landmarks,identity=M.identity(),rows=[];
 for(const [name,parent,child,pivot] of [['SC','thorax','clavicle',L.SC],['AC','clavicle','scapula',L.AC],['GH','scapula','humerus',L.GH],['elbow','humerus','ulna',L.EL],['PRUJ-axis','ulna','radius',radialHead],['DRUJ-axis','ulna','radius',distalUlna]]){
  if(!pivot||!bones[child])continue;
  const a=bones[parent]||identity,b=bones[child],pa=apply(a,pivot),pb=apply(b,pivot);
  rows.push({joint:name,parent,child,restPivot:pivot,parentMeshLocal:sub(pivot,a.rest),childMeshLocal:sub(pivot,b.rest),parentWorld:pa,childWorld:pb,sharedPivotErrorMm:M.len(sub(pa,pb))/model.scale,articularSurfaceGapMm:null});
 }
 return rows;
}
function shoulderMotion(model,pose,input){
 const degrees=x=>x*180/Math.PI,relative=q=>qmul(F.conjugate(model.thoraxFrame),q),st=q=>F.yxz(relative(qmul(q,model.scapulaFrame)));
 const current=st(pose.qS),rest=st([0,0,0,1]),prior=model.pose({...input,fit:false,backend:'reference'}),reference=st(prior.qS);
 const named=a=>({internalRotation:degrees(a[0]),upwardRotation:-degrees(a[1]),posteriorTilt:degrees(a[2])});
 const delta=(a,b)=>a.map((x,i)=>Math.atan2(Math.sin(x-b[i]),Math.cos(x-b[i])));
 const shaft=norm(sub(model.landmarks.EL,model.landmarks.GH)),world=qrot(pose.qH,shaft),thoraxUp=qrot(model.thoraxFrame,[0,1,0]);
 const scapulaUp=qrot(qmul(pose.qS,model.scapulaFrame),[0,1,0]);
 return {coordinateSystem:'canonical-right anatomical thorax; ST Y-X-Z decomposition',units:'degrees',
  humerothoracicElevation:degrees(Math.acos(M.clamp(-dot(world,thoraxUp),-1,1))),
  scapulothoracic:named(current),scapulothoracicChangeFromRest:named(delta(current,rest)),
  unconstrainedPriorScapulothoracic:named(reference),contactChangeFromPrior:named(delta(current,reference)),
  humeralProximalAxisToScapularYAxis:degrees(Math.acos(M.clamp(-dot(world,scapulaUp),-1,1))),
  contributionsAreNotAdditive:true,clinicalROMorDiseaseClassification:null};
}
function muscle(state){
 const B=state.binding,length=c=>{let n=0,p=root.FiberPath.curve(c,0);for(let i=1;i<=32;i++){const q=root.FiberPath.curve(c,i/32);n+=M.len(sub(q,p));p=q}return n/B.definition.source.worldUnitsPerMm};
 return {id:B.definition.id,curves:B.field.controls.map((c,i)=>{const restMm=length(c),currentMm=length(state.posed[i]);return {lane:i,referenceLengthMm:restMm,currentLengthMm:currentMm,ratio:currentMm/restMm}}),attachments:B.definition.attachments.map(a=>({bone:a.bone,vertices:a.indices.length,confidence:a.confidence,coordinateFrame:'mesh-rest-atlas',units:'world'})),activation:null,force:null,anatomicalViaPointsVerified:false};
}
root.RigDiagnostics=Object.freeze({longBone,matrix,frames,connections,muscle,shoulderMotion});
})(globalThis);
