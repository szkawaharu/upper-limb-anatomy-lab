import numpy as np
from constrained_cage import ConstrainedCage
from conforming_surface import conform
p=np.array([[0.,0,0],[1,0,0],[0,1,0],[.2,.2,1],[.2,.2,-1]])
t=np.array([[0,1,2,3],[0,2,1,4]]);classes=np.array([0,0,0,-1,1]);c=ConstrainedCage(p,t,classes)
rng=np.random.default_rng(72);b=rng.dirichlet(np.ones(4),100);points=np.r_[b@p[t[0]],b@p[t[1]],p,np.array([[.2,.2,0]])]
ids,weights,which=c.embed(points);assert np.max(np.abs(c.positions(p,ids,weights)-points))<1e-12
# The internal fixed facet is emitted once, has no free supports, and stays
# exactly attached while the adjacent volume control is free to move.
r=conform(c,p,np.array([[0,1,2]]),np.array([0,0,0,.5,1]));assert len(r['faces'])==1
cl,ids,bary,_=c.fixed_patches(r['points'],r['weights']);moved=p.copy();moved[3]+=[.1,.2,.3];now=c.positions(moved,ids,bary);assert np.max(np.abs(now[:3]-p[:3]))<1e-12
assert np.linalg.norm(now[3]-p[3])>.1
print('PASS constrained cage: independent two-cell barycentric reconstruction, exact internal fixed area, no duplicate coplanar facet, adjacent free tissue remains movable')
