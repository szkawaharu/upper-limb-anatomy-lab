import assert from 'node:assert/strict';
import '../engine/math.js';
import '../engine/deformation-audit.js';
import '../scripts/experimental/deformation-graph.js';
const M=AnatomyMath,D=DeformationGraph;
// Known rigid covariance (including almost planar data) has a known answer;
// compare transformed basis vectors, never quaternion sign conventions.
for(const r of [0,.7,Math.PI-1e-6,Math.PI,Math.PI+1e-6,5.8]){
 const q=M.qaxis([.3,.4,.5],r),H=new Float64Array(9);
 for(const p of [[1,0,0],[0,1,0],[0,0,.0001]]){const x=M.qrot(q,p);for(let i=0;i<3;i++)for(let j=0;j<3;j++)H[i*3+j]+=p[i]*x[j];}
 const R=D.rotation(H);
 for(let i=0;i<3;i++)for(let j=0;j<3;j++){const p=[0,0,0];p[j]=1;assert(Math.abs(R[i*3+j]-M.qrot(q,p)[i])<1e-10);}
}
const rest=Float32Array.of(0,0,0,1,0,0,0,1,0,0,0,1),faces=Uint32Array.of(0,2,1,0,1,3,0,3,2,1,2,3),edges=[0,1,0,2,0,3,1,2,1,3,2,3],edgeWeights=Array(6).fill(1),embed=new Float32Array(32);
for(let i=0;i<4;i++){embed[i*8]=i;embed[i*8+4]=1;}
const G={nodes:[...rest],nodeWeights:[0,1,.4,.6],nodeClass:[0,1,2,2],edges,edgeWeights,free:[2,3],fixed:[0,1],cellWorld:1,factor:{lowerOffsets:[0,0,1],lowerIndices:[0],lowerValues:[-1/Math.sqrt(3)],upperOffsets:[0,1,1],upperIndices:[1],upperValues:[-1/Math.sqrt(3)],diagonal:[Math.sqrt(3),Math.sqrt(8/3)]}};
const B={rest,faces,weights:Float32Array.of(0,1,.4,.6),graph:D.prepare(G,embed),definition:{id:'closed-tetrahedron',source:{worldUnitsPerMm:1}}};
// A double-precision Cholesky/eigensolve has round-off even for zero motion.
const restError=p=>Math.max(...p.map((x,i)=>Math.abs(x-rest[i])));
const identity=M.identity(),neutral=D.pose(B,identity,identity),neutralPoints=D.deform(neutral);assert(restError(neutralPoints)<1e-12);
const rigid={q:M.qaxis([.3,.4,.5],2.1),rest:[0,0,0],p:[.4,.5,.6]},rigidState=D.pose(B,rigid,rigid),rigidPoints=D.deform(rigidState);
for(let i=0;i<4;i++)assert(M.len(M.sub(rigidPoints.subarray(i*3,i*3+3),M.apply(rigid,rest.subarray(i*3,i*3+3))))<1e-7);
assert.deepEqual([...D.deform(neutral)],[...neutralPoints],'Cached other-arm pose changed after a new solve');
B.graphQuadrature={sourceVertices:[0,1,2,3],faces:[...faces],multipliers:[1,1,1,1],controlDet:[0,0,0,1],referenceVolume:1/6,sourceCenter:[0,0,0],pureBVector:[0,0,0]};
const bent={q:M.qaxis([0,0,1],.9),p:[.8,.2,.1],rest:[1,0,0]},state=D.pose(B,identity,bent,{volume:true}),now=D.deform(state),report=DeformationAudit.audit(state,now);
assert(report.topology.closed);assert.equal(report.nonfinite,0);assert.equal(report.stretch,null);assert.equal(report.localDeterminantApprox,null);assert(report.attachmentMaxMm<1e-7);assert(Math.abs(report.volumeRatio-1)<1e-5);
assert(M.len(M.sub(D.point(state,0),M.apply(identity,[0,0,0])))<1e-12);assert(M.len(M.sub(D.point(state,1),M.apply(bent,[1,0,0])))<1e-12);
// Independent finite differences of the closed integral test the constraint
// derivative, including the attachment-boundary basis, rather than the solve.
for(const boundary of [null,'quartic-dirichlet-bubble']){
 B.graphBoundary=boundary;const frame=D.pose(B,identity,bent),gradient=D.volumeData(frame,true).grad;
 for(const i of frame.graph.free)for(let k=0;k<3;k++){
  const j=i*3+k,v=frame.graph.x[j],h=1e-6;frame.graph.x[j]=v+h;const hi=D.volumeData(frame).volume;frame.graph.x[j]=v-h;const lo=D.volumeData(frame).volume;frame.graph.x[j]=v;
  assert(Math.abs((hi-lo)/(2*h)-gradient[j])<1e-9,'Source volume derivative disagrees with finite differences');
 }
}
B.graphBoundary=null;
const returned=D.pose(B,identity,identity);assert(restError(D.deform(returned))<1e-12);
console.log('Experimental ARAP: known rigid and planar rotations, fixed boundaries, closed-volume projection, stateless return, cached bilateral isolation. No source-muscle shape/medical certification.');
