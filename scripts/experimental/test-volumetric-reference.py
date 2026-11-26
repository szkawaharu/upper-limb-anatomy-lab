"""Independent continuum checks on a closed unit tetrahedron."""
import numpy as np
from volumetric_reference import Solid
p=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1.]],float)
solid=Solid(p,[[0,1,2,3]],[0])
q=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
E,g,J=solid.evaluate(p@q.T+[.2,-.3,.4])
assert abs(E)<1e-12 and np.max(np.abs(g))<1e-12 and np.max(np.abs(J-1))<1e-12
x=p@np.array([[1.1,.2,.1],[.1,.9,.2],[0,.1,1.05]]).T
E,g,J=solid.evaluate(x)
for i in range(4):
 for k in range(3):
  hi=x.copy();lo=x.copy();hi[i,k]+=1e-6;lo[i,k]-=1e-6
  derivative=(solid.evaluate(hi,gradient=False)[0]-solid.evaluate(lo,gradient=False)[0])/2e-6
  assert abs(derivative-g[i,k])<1e-7,(i,k,derivative,g[i,k])
bad=p.copy();bad[3,2]=-1
assert not np.isfinite(solid.evaluate(bad)[0])
scaled=p*2;assert abs(solid.evaluate(scaled)[2][0]-8)<1e-12
start=p.copy();start[1:]*=1.1
result,r=solid.minimize(start,iterations=300,residual_tolerance_mm=1e-7)
assert r['converged'] and r['nonPositiveElements']==0 and abs(r['volumeRatio']-1)<1e-6,r
for bulk in [3,3000]:
 result,r=solid.minimize(start,bulk=bulk,iterations=1000,residual_tolerance_mm=1e-7)
 assert r['converged'] and r['preconditionedGradientMaxMm']<1e-7,r
print('PASS continuum: rigid invariance, analytic stress finite differences, orientation rejection, 8x volume, constrained relaxation; no source-muscle validation implied')
