import numpy as np
from attachment_refinement import refine
from constrained_cage import ConstrainedCage
from volumetric_reference import Solid
p=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1]])
t=np.array([[0,1,2,3]]);cl=np.zeros(4,int);declared=np.array([[0,1,2]]);edges=np.array([[0,1],[1,2],[2,0]])
n,c,labels,r,parent=refine(p,t,cl,declared,edges);assert np.array_equal(n[:4],p)and np.all(labels[:4]==0)and np.all(parent==0)
keys=set(tuple(sorted(f))for a in c for f in [a[[0,1,2]],a[[0,1,3]],a[[0,2,3]],a[[1,2,3]]]);assert (0,1,2)in keys
C=ConstrainedCage(n,c,labels);q=np.array([[.2,.2,0],[.2,0,.2],[0,0,.5],[.2,.2,.2]]);ids,bary,_=C.embed(q);now=n.copy();now[labels<0]+=[.01,.02,.03];result=C.positions(now,ids,bary)
assert np.max(np.abs(result[0]-q[0]))<1e-12
assert np.all(np.linalg.norm(result[1:]-q[1:],axis=1)>.005)
volume=np.linalg.det(np.stack([n[c[:,k]]-n[c[:,0]]for k in [1,2,3]],axis=2))/6
old=Solid(p,t,[0,1,2,3]);new=Solid(n,c,np.flatnonzero(labels==0),element_weights=volume/volume.sum());A=np.array([[1.1,.2,0],[0,.85,.03],[0,0,1.05]])
assert abs(old.evaluate(p@A.T,30)[0]-new.evaluate(n@A.T,30)[0])<1e-12
print('PASS attachment refinement: exact declared fixed facet, artificial fixed face/edge/interior freed, original controls unchanged, total volume and positive reference cells retained')
