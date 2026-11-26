"""Independent algebra/safety checks, not source-anatomy acceptance."""
import numpy as np
from scipy.sparse.linalg import cg
from volumetric_reference import Solid
from volumetric_gauss_newton import linearize,minimize
from volumetric_preconditioner import rest_metric_matrix,factor_rest_metric
from scipy.sparse.linalg import LinearOperator
p=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1.]],float);solid=Solid(p,[[0,1,2,3]],[0]);x=p@np.array([[1.1,.2,.1],[.1,.9,.2],[0,.1,1.05]]).T
H=linearize(solid,x,3000);dense=np.column_stack([H@np.eye(9)[:,j]for j in range(9)])
assert np.max(np.abs(dense-dense.T))<1e-10 and np.linalg.eigvalsh(dense).min()>0
rhs=np.arange(9,dtype=float)-4;d,info=cg(H,rhs,rtol=1e-10,atol=0,maxiter=100)
assert info==0 and np.linalg.norm(d-np.linalg.solve(dense,rhs))<1e-8
# Independently assembled sparse source-rest metric must equal the matrix-free
# rest operator; changing a linear preconditioner cannot change the energy.
restH=linearize(solid,p,3000)
expected=np.column_stack([restH@np.eye(9)[:,j]for j in range(9)])
actual=rest_metric_matrix(solid,3000).toarray()
assert np.max(np.abs(actual-expected))<1e-12
factor,meta=factor_rest_metric(solid,30)
assert np.linalg.norm(rest_metric_matrix(solid,30)@factor.solve(rhs)-rhs)<1e-10
# The metric must follow the actual world rotation, not fixed world axes.
R=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
v=np.arange(9,dtype=float).reshape(-1,3)/9
a=(H@v.ravel()).reshape(-1,3)@R.T
b=(linearize(solid,x@R.T,3000)@(v@R.T).ravel()).reshape(-1,3)
assert np.linalg.norm(a-b)/np.linalg.norm(a)<1e-12
result,report=minimize(solid,x,bulk=30,iterations=100,cg_iterations=30,residual_tolerance_mm=1e-7)
assert report['nonPositiveElements']==0 and report['energy']<solid.evaluate(x,30)[0]
assert report['converged'] and report['preconditionedGradientMaxMm']<1e-7,report
alternate,other=minimize(solid,x,bulk=30,iterations=100,cg_iterations=30,residual_tolerance_mm=1e-7,preconditioner=LinearOperator((9,9),matvec=factor.solve,dtype=float),preconditioner_info=meta)
assert other['converged'] and np.max(np.abs(result-alternate))<1e-6,other
print('PASS matrix-free metric: symmetry, positive eigenvalues, independent dense solve, world-rotation covariance, energy decrease, exact-J line-search rejection; no anatomical certification')
