import numpy as np
from point_cage import PointCage
from conforming_surface import conform
p=np.array([[.17,.13,.23],[3.37,.22,.41],[.32,2.91,.49],[.22,.34,3.15]])
w=np.array([0,0,.4,1]);c=PointCage(p,w,.8);classes,ids,bary,tet=c.fixed_patches(p,w);assert (classes>=0).sum()==3
assert np.max(np.abs(c.positions(c.nodes,ids,bary)-p))<1e-12
# Move all free cage controls. Only the three declared fixed points stay put;
# unlike grid-support pinning, the neighboring muscle point remains movable.
now=c.nodes.copy();now[classes<0]+=[.01,.02,.03];y=c.positions(now,ids,bary)
assert np.max(np.abs(y[np.isin(w,[0,1])]-p[np.isin(w,[0,1])]))<1e-12
assert np.linalg.norm(y[2]-p[2])>.01
f=np.array([[0,2,1],[0,1,3],[1,2,3],[2,0,3]])
r=conform(c,p,f,w);assert r['sourceAreaMaxRelativeError']<1e-9
for _,face,b in r['newVertexProvenance']:assert np.isclose(sum(b),1)
print('PASS exact attachment control cage: no neighboring-support pinning, source reconstruction, affine surface subdivision including boundary-aligned facets')
