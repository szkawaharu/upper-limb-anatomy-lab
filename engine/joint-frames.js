/* Explicit rigid transforms: mesh coordinates, anatomical frames and joint
 * pivots are separate. Rotations are proper quaternions, never reflections. */
(function(root){
'use strict';
const M=root.AnatomyMath,{add,sub,mul,dot,cross,norm,qmul,qrot,apply}=M;
const conjugate=q=>[-q[0],-q[1],-q[2],q[3]];
function basis(x,y,z){
 if(![...x,...y,...z].every(Number.isFinite)||Math.abs(dot(x,cross(y,z))-1)>1e-6||[x,y,z].some(v=>Math.abs(dot(v,v)-1)>1e-6)||Math.max(Math.abs(dot(x,y)),Math.abs(dot(y,z)),Math.abs(dot(x,z)))>1e-6)throw Error('Degenerate anatomical frame');
 const R=[[x[0],y[0],z[0]],[x[1],y[1],z[1]],[x[2],y[2],z[2]]],t=R[0][0]+R[1][1]+R[2][2];let q;
 if(t>0){const s=Math.sqrt(1+t)*2;q=[(R[2][1]-R[1][2])/s,(R[0][2]-R[2][0])/s,(R[1][0]-R[0][1])/s,s/4]}
 else{let i=0;for(let j=1;j<3;j++)if(R[j][j]>R[i][i])i=j;const j=(i+1)%3,k=(i+2)%3,s=Math.sqrt(1+R[i][i]-R[j][j]-R[k][k])*2;q=[0,0,0,0];q[i]=s/4;q[j]=(R[j][i]+R[i][j])/s;q[k]=(R[k][i]+R[i][k])/s;q[3]=(R[k][j]-R[j][k])/s}
 const n=Math.hypot(...q);return q.map(x=>x/n);
}
function yxz(q){
 const x=qrot(q,[1,0,0]),y=qrot(q,[0,1,0]),z=qrot(q,[0,0,1]);
 return [Math.atan2(z[0],z[2]),Math.asin(M.clamp(-z[1],-1,1)),Math.atan2(x[1],y[1])];
}
function fromYXZ(a){return qmul(M.qaxis([0,1,0],a[0]),qmul(M.qaxis([1,0,0],a[1]),M.qaxis([0,0,1],a[2])))}
function clavicle(SC,AC,up=[0,1,0]){const z=norm(sub(AC,SC)),x=norm(cross(up,z)),y=norm(cross(z,x));return basis(x,y,z)}
function thorax(IJ,C7,PX,T8){const upper=mul(add(IJ,C7),.5),lower=mul(add(PX,T8),.5),y=norm(sub(upper,lower));let z=norm(cross(sub(C7,IJ),sub(lower,IJ)));if(dot(z,[1,0,0])<0)z=mul(z,-1);const x=norm(cross(y,z));return basis(x,norm(cross(z,x)),z)}
function scapula(TS,AA,AI){const z=norm(sub(AA,TS));let x=norm(cross(sub(TS,AA),sub(AI,AA)));if(dot(x,[0,0,-1])<0)x=mul(x,-1);return basis(x,norm(cross(z,x)),z)}
function swing(a,b){const k=cross(a,b),d=dot(a,b);if(d<-.999999)throw Error('Ambiguous neutral swing');const q=[...k,1+d],n=Math.hypot(...q);return q.map(x=>x/n)}
// A child is attached at a shared world-rest point, not at its mesh origin.
function attach(parent,pivot,relativeQ,renderOrigin=pivot){
 const q=qmul(parent.q,relativeQ),joint=apply(parent,pivot);
 return {q,p:add(joint,qrot(q,sub(renderOrigin,pivot))),rest:renderOrigin};
}
root.JointFrames=Object.freeze({conjugate,basis,yxz,fromYXZ,clavicle,scapula,thorax,swing,attach});
})(globalThis);
