import numpy as np
from volumetric_reference import Solid
from projected_solid_metric import assemble
from cage_solver import minimize
p=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1],[.25,.25,.25]])
t=np.array([[4,1,2,3],[0,4,2,3],[0,1,4,3],[0,1,2,4]]);s=Solid(p,t,[0,1,2,3]);F=np.array([[1.1,.1,0],[0,.85,.03],[0,0,1.05]])
x=p@F.T;x[4]+=[.025,-.015,.01];H,_=assemble(s,x,30,False);fd=[]
for k in range(3):
 hi=x.copy();lo=x.copy();hi[4,k]+=1e-6;lo[4,k]-=1e-6;fd.append(((s.evaluate(hi,30)[1]-s.evaluate(lo,30)[1])/2e-6)[4])
assert np.max(np.abs(H.toarray()-np.array(fd).T))<1e-5
now,r=minimize(s,x,tolerance_mm=1e-9);assert r['converged']and np.max(np.abs(now-p@F.T))<1e-8
assert np.array_equal(now[s.fixed],x[s.fixed]);assert all(a['minJ']>0 for a in r['rows'])
assert all(a['energy']>=b['energy']for a,b in zip(r['rows'],r['rows'][1:]))
print('PASS safeguarded global Newton: independent non-rest stress derivative, analytic homogeneous-deformation equilibrium, fixed boundary, monotone energy, positive J')
