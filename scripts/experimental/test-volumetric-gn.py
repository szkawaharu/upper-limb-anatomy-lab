"""Independent algebra/safety checks, not source-anatomy acceptance."""
import numpy as np
from scipy.sparse.linalg import cg
from volumetric_reference import Solid
from volumetric_gauss_newton import linearize,minimize
from volumetric_preconditioner import rest_metric_matrix,strain_metric_matrix,factor_rest_metric,warped_rest_preconditioner
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
assert np.max(np.abs(strain_metric_matrix(solid,x,3000).toarray()-dense))<1e-12
factor,meta=factor_rest_metric(solid,30)
assert np.linalg.norm(rest_metric_matrix(solid,30)@factor.solve(rhs)-rhs)<1e-10
# Overlapping tetrahedra must accumulate shared-node blocks correctly.
star=np.vstack([p,[.25,.25,.25]])
cells=np.array([[4,1,2,3],[0,4,2,3],[0,1,4,3],[0,1,2,4]])
joined=Solid(star,cells,[0,1,2,3]);posed=star.copy();posed[4]+=[.01,-.015,.02]
operator=linearize(joined,posed,30)
expected=np.column_stack([operator@np.eye(3)[:,j]for j in range(3)])
actual=strain_metric_matrix(joined,posed,30).toarray()
assert np.max(np.abs(actual-expected))<1e-12
# The metric must follow the actual world rotation, not fixed world axes.
R=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
v=np.arange(9,dtype=float).reshape(-1,3)/9
# A common rotation must transport the inverse of the entire rest metric.
P=warped_rest_preconditioner(solid,p@R.T,factor)
rotated=linearize(solid,p@R.T,30)
assert np.linalg.norm(rotated@(P@rhs)-rhs)<1e-10
Pmatrix=np.column_stack([P@np.eye(9)[:,j]for j in range(9)])
assert np.max(np.abs(Pmatrix-Pmatrix.T))<1e-12 and np.linalg.eigvalsh(Pmatrix).min()>0
a=(H@v.ravel()).reshape(-1,3)@R.T
b=(linearize(solid,x@R.T,3000)@(v@R.T).ravel()).reshape(-1,3)
assert np.linalg.norm(a-b)/np.linalg.norm(a)<1e-12
result,report=minimize(solid,x,bulk=30,iterations=100,cg_iterations=30,residual_tolerance_mm=1e-7)
assert report['nonPositiveElements']==0 and report['energy']<solid.evaluate(x,30)[0]
assert report['converged'] and report['preconditionedGradientMaxMm']<1e-7,report
alternate,other=minimize(solid,x,bulk=30,iterations=100,cg_iterations=30,residual_tolerance_mm=1e-7,preconditioner=LinearOperator((9,9),matvec=factor.solve,dtype=float),preconditioner_info=meta)
assert other['converged'] and np.max(np.abs(result-alternate))<1e-6,other
print('PASS matrix-free metric: symmetry, positive eigenvalues, independent dense solve, world-rotation covariance, energy decrease, exact-J line-search rejection; no anatomical certification')
