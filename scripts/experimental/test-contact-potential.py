import numpy as np
from contact_potential import PlaneContact,AugmentedSolid
from volumetric_reference import Solid
from projected_solid_metric import assemble
p=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1],[.25,.25,.25]]);t=np.array([[4,1,2,3],[0,4,2,3],[0,1,4,3],[0,1,2,4]]);s=Solid(p,t,[0,1,2,3])
vertices=np.array([[0,1,2,4],[1,2,3,4]]);bary=np.array([[.1,.1,.1,.7],[.1,.2,.1,.6]]);n=np.array([[1.,0,0],[0,1,0]]);offsets=np.array([.6,.7]);w=np.array([.3,.7]);potential=PlaneContact(5,vertices,bary,n,offsets,w,20);aug=AugmentedSolid(s,[potential]);energy,g,_=aug.evaluate(p);H,_=assemble(aug,p,30,False)
fd=[]
for k in range(3):
 hi=p.copy();lo=p.copy();hi[4,k]+=1e-6;lo[4,k]-=1e-6
 assert abs((aug.evaluate(hi,gradient=False)[0]-aug.evaluate(lo,gradient=False)[0])/2e-6-g[4,k])<1e-7
 fd.append(((aug.evaluate(hi)[1]-aug.evaluate(lo)[1])/2e-6)[4])
assert np.max(np.abs(H.toarray()-np.array(fd).T))<1e-5
q=p.copy();q[:]+=[2,2,0];assert potential.evaluate(q)[0]==0 and np.max(np.abs(potential.evaluate(q)[1]))==0
R=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
rot=PlaneContact(5,vertices,bary,n@R.T,offsets,w,20);assert abs(rot.evaluate(p@R.T)[0]-potential.evaluate(p)[0])<1e-12
assert np.max(np.abs(rot.evaluate(p@R.T)[1]-potential.evaluate(p)[1]@R.T))<1e-12
print('PASS contact potential: independent energy/stress/Hessian derivatives, unilateral inactive zero, material-point quadrature, rigid-rotation covariance')
