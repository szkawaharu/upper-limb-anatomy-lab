"""Independent rest derivative and full/reduced projection checks."""
import numpy as np
from volumetric_reference import Solid
from linear_elastic_modes import RestElasticity,build_modes
from vector_reduced_solid import VectorReducedSolid
from volumetric_gauss_newton import linearize
p=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[.25,.25,.25],[9,9,9]],float)
t=np.array([[4,1,2,3],[0,4,2,3],[0,1,4,3],[0,1,2,4]])
assert np.all(np.linalg.det(np.stack([p[t[:,i]]-p[t[:,0]]for i in [1,2,3]],axis=2))>0)
s=Solid(p,t,[0,1,2,3]);elastic=RestElasticity(s,30)
v=np.arange(18,dtype=float).reshape(-1,3)/30
for bulk in [3,30,3000]:
 H=RestElasticity(s,bulk);eps=1e-6
 fd=(s.evaluate(p+eps*v,bulk)[1]-s.evaluate(p-eps*v,bulk)[1])/(2*eps)
 assert np.linalg.norm(H.apply(v)-fd)/np.linalg.norm(fd)<1e-8
# Infinitesimal common rigid motion creates no rest strain or force.
W=np.array([[0,-.3,.2],[.3,0,-.1],[-.2,.1,0]])
assert np.linalg.norm(elastic.apply(p@W.T+np.array([.2,-.1,.4])))<1e-12
boundary=np.zeros((4,3));boundary[3]=[.02,-.01,.03]
x,report=elastic.solve_boundary(boundary,rtol=1e-12)
assert report['cgInfo']==0 and report['preconditionedResidualMm']<1e-12
assert np.array_equal(x[s.fixed],boundary)
# Compare with a separately assembled dense exact rest tangent.
D=np.column_stack([elastic.operator@np.eye(3)[:,k]for k in range(3)])
initial=np.zeros_like(p);initial[s.fixed]=boundary
expected=np.linalg.solve(D,-elastic.apply(initial)[s.free].ravel())
assert np.max(np.abs(expected-x[s.free].ravel()))<1e-12
lifting=np.zeros((len(p),4));lifting[3]=np.r_[p[3],1]
lifting[s.free]=-s.factor.solve(s.coupling@lifting[s.fixed])
r=build_modes(s,lifting,np.array([0,0,0,1]),30);V=r['vectorBasis']
assert V.shape==(6,3,3)and np.all(V[s.fixed]==0)and np.all(V[5]==0)
model=VectorReducedSolid(p,t,s.fixed,V);base=p.copy();base[3]+=[.02,-.01,.03];base=s.predict_boundary(p,base[s.fixed]);q=np.arange(3)*.0002
E,g,J,H=model.evaluate(base,q,30,hessian=True);X=model.positions(base,q)
assert abs(E-s.evaluate(X,30)[0])<1e-12
for k in range(len(q)):
 a=q.copy();b=q.copy();a[k]+=1e-6;b[k]-=1e-6
 fd=(model.evaluate(base,a,30)[0]-model.evaluate(base,b,30)[0])/2e-6
 assert abs(g[k]-fd)<1e-8
fullH=linearize(s,X,30);Rfree=V[s.free].reshape(-1,V.shape[2]);expected=Rfree.T@np.column_stack([fullH@Rfree[:,k]for k in range(V.shape[2])])
assert np.linalg.norm(H-expected)/np.linalg.norm(expected)<1e-12
X,q,report=model.minimize(base,iterations=100,tolerance_mm=1e-8)
assert report['converged']and report['fullReferenceConverged'],report
assert np.array_equal(X[s.fixed],base[s.fixed])
print('PASS exact rest elasticity: independent stress derivative, infinitesimal rigid null mode, boundary response/dense solve, rank, zero attachment modes; vector reduction: independent energy gradient/full metric and complete-space equilibrium')
