import numpy as np
from scipy.spatial.transform import Rotation
from reduced_volume_basis import build_basis
from reduced_volume_model import ReducedSolid,affine_boundary_lift
from volumetric_gauss_newton import linearize
p=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1.],[4,4,4]],float);t=np.array([[0,1,2,3]]);w=np.array([0,1,.3,.7]);b=build_basis(p,t,w,controls=2);phi=b['basis'];model=ReducedSolid(p,t,b['fixed'],phi)
assert b['discardedNumericalModes']==6 and phi.shape==(5,2)
assert np.all(phi[b['fixed']]==0)and np.all(phi[4]==0)and model.metricError<1e-12
assert np.array_equal(b['lifting'][0],np.zeros(4))
assert np.array_equal(b['lifting'][1],np.r_[p[1],1])
# Affine lifting is exact on both fixed regions and under common rigid motion.
R=Rotation.from_rotvec([.2,.3,.4]).as_matrix();at=np.array([.3,-.4,.2]);base=p@R.T+at
assert np.array_equal(affine_boundary_lift(p,b['lifting'],R,at,R,at),base)
# With an identity reference bone, the four-channel lift equals an
# independently solved full-FEM harmonic displacement predictor.
RB=Rotation.from_rotvec([0,.02,.03]).as_matrix();bt=np.array([.01,.02,.03]);lift=affine_boundary_lift(p,b['lifting'],np.eye(3),np.zeros(3),RB,bt)
target=p[b['fixed']].copy();target[1]=RB@p[1]+bt
predicted=model.solid.predict_boundary(p,target)
assert np.max(np.abs(lift[:4]-predicted[:4]))<1e-12
E,g,J,H=model.evaluate(base,np.zeros((2,3)),30,hessian=True)
assert abs(E)<1e-12 and np.linalg.norm(g)<1e-12
base=p.copy();base[1]=[1,.1,.1];q=np.arange(6,dtype=float).reshape(2,3)*.002
E,g,J,H=model.evaluate(base,q,30,hessian=True)
for i in range(2):
 for j in range(3):
  hi=q.copy();lo=q.copy();hi[i,j]+=1e-6;lo[i,j]-=1e-6
  fd=(model.evaluate(base,hi,30)[0]-model.evaluate(base,lo,30)[0])/2e-6
  assert abs(fd-g[i,j])<1e-7
# Independently project the full matrix-free metric onto the reduced basis.
x=base+phi@q;fullH=linearize(model.solid,x,30);Rfree=np.kron(phi[model.solid.free],np.eye(3));expected=Rfree.T@np.column_stack([fullH@Rfree[:,j]for j in range(Rfree.shape[1])])
assert np.linalg.norm(expected-H)/np.linalg.norm(H)<1e-12
x,q,r=model.minimize(base,bulk=30,iterations=100,tolerance_mm=1e-7)
assert r['converged']and r['fullReferenceConverged']and r['nonPositiveElements']==0,r
assert np.array_equal(x[b['fixed']],base[b['fixed']])
print('PASS reduced solid: exact pins, inactive node, rank/whitening, rigid invariance, finite-difference gradient, independent full-metric projection, complete-space continuum convergence')
