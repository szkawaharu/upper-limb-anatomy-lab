import numpy as np
from tetrahedral_cage import GridCage
p=np.random.default_rng(42).uniform(-3,4,(100,3));c=GridCage(p,1.3);ids,bary,t=c.embed(p)
assert np.max(np.abs(c.positions(c.nodes,ids,bary)-p))<1e-12
assert bary.min()>-1e-12 and np.max(np.abs(bary.sum(1)-1))<1e-12
A=np.array([[1.2,.2,.1],[-.1,.9,.2],[0,.3,1.1]]);b=np.array([3,-2,4])
assert np.max(np.abs(c.positions(c.nodes@A.T+b,ids,bary)-(p@A.T+b)))<1e-12
q=c.nodes[c.tets];E=np.stack([q[:,k]-q[:,0]for k in [1,2,3]],axis=2);assert np.all(np.linalg.det(E)>0)
# Grid-boundary tie handling and opposite original footprint patches.
points=np.array([[0,0,0],[0,.2,.3],[5,0,0],[5,.2,.3]])
c=GridCage(points,1);classes,ids,bary,t=c.fixed_patches(points,[0,0,1,1]);nodes=c.nodes.copy();nodes[classes==0]+=[0,1,2];nodes[classes==1]+=[2,1,0];x=c.positions(nodes,ids,bary)
assert np.max(np.abs(x-points-np.array([[0,1,2],[0,1,2],[2,1,0],[2,1,0]])))<1e-12
print('PASS cage: source reconstruction, affine reproduction, orientation, grid ties, exact wide source footprint supports')
