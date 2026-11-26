import assert from'node:assert/strict';
import '../engine/math.js';import '../engine/muscle-deformation.js';import '../engine/deformation-audit.js';
const rest=Float32Array.of(0,0,0,.001,0,0,0,.001,0,0,0,.001),faces=Uint32Array.of(0,2,1,0,1,3,0,3,2,1,2,3),B={rest,faces,weights:new Float32Array(4),lane:new Float32Array(4),material:new Float32Array(4),definition:{id:'unit-tetrahedron',source:{worldUnitsPerMm:.001}}},state={binding:B,a:AnatomyMath.identity(),b:AnatomyMath.identity()},now=Float32Array.from(rest,x=>2*x),r=DeformationAudit.audit(state,now);
assert(r.topology.closed);assert(Math.abs(r.volumeRatio-8)<1e-10);assert.equal(r.edgeRatio[1],2);assert.equal(r.areaRatio[1],4);assert.equal(r.worstEdge.vertices.length,2);assert.equal(r.worstArea.vertices.length,3);assert.equal(r.worstEdge.ratio,2);assert.equal(r.worstArea.ratio,4);
const [a,b]=r.worstEdge.vertices;assert(Math.abs(r.worstEdge.restLengthMm-Math.hypot(...rest.slice(a*3,a*3+3).map((x,i)=>x-rest[b*3+i]))/.001)<1e-12);assert(r.worstArea.restAreaMm2>=.5);
console.log('Independent scaled closed tetrahedron: physical edge/area units, outlier identities, volume ratio');
