import numpy as np
from volumetric_reference import Solid
from projected_solid_metric import assemble
p=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1.],[.25,.25,.25]])
t=np.array([[4,1,2,3],[0,4,2,3],[0,1,4,3],[0,1,2,4]]);s=Solid(p,t,[0,1,2,3]);H,info=assemble(s,p,30);H=H.toarray();fd=[]
for k in range(3):
 hi=p.copy();lo=p.copy();hi[4,k]+=1e-6;lo[4,k]-=1e-6;fd.append(((s.evaluate(hi,30)[1]-s.evaluate(lo,30)[1])/2e-6)[4])
assert np.max(np.abs(H-np.array(fd).T))<1e-6
assert np.linalg.eigvalsh(H).min()>0
R=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
rotated,_=assemble(s,p@R.T,30);assert np.max(np.abs(rotated.toarray()-R@H@R.T))<1e-10
print('PASS projected sparse tangent: independent shared-element stress derivative at rest, positivity, rigid-rotation covariance')
