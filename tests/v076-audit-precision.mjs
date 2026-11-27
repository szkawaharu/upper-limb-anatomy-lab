import assert from'node:assert/strict';
await import('../engine/math.js');await import('../engine/deformation-audit.js');
const A=globalThis.DeformationAudit,faces=[0,2,1,0,1,3,0,3,2,1,2,3],tet=[[0,0,0],[1,0,0],[0,1,0],[0,0,1]];
for(const scale of [1e-12,1e-8,1e-4,1,1e4,1e8]){
 const points=tet.map(p=>p.map(x=>x*scale)),topology=A.topology(points.flat(),faces);assert(topology.closed);assert.equal(topology.degenerate,0);
 assert.equal(topology.smallAbsoluteAreaFaces,scale<=1e-8?4:0);
 assert(A.triangleQuality(...points.slice(0,3)).reliablyNoncollinear);
 const line=[[0,0,0],[1,1,1],[2,2,2]].map(p=>p.map(x=>x*scale));assert(!A.triangleQuality(...line).reliablyNoncollinear);
}
const tiny=[[1,2,3],[1+1e-7,2,3],[1,2+1e-7,3]].map(p=>Array.from(Float32Array.from(p)));
assert(!A.triangleQuality(...tiny).reliablyNoncollinear,'Actual Float32 coordinate collapse must still be rejected');
const valid=[[1,2,3],[1+1e-6,2,3],[1,2+1e-6,3]].map(p=>Array.from(Float32Array.from(p)));
assert(A.triangleQuality(...valid).reliablyNoncollinear);
assert(!A.triangleQuality([NaN,0,0],[1,0,0],[0,1,0]).reliablyNoncollinear);
console.log('PASS scale-independent mesh rank: small valid tetrahedra survive unit conversion; collinear, collapsed Float32 and nonfinite faces remain rejected; old small-area counts retained');
