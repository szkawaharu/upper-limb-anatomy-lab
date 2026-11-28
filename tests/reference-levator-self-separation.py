"""Synthetic check of finite-pair separation, rigid transport and derivatives."""
import sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experimental'))
from local_self_separation import SelfSeparation,reference_axis,transported_axis
from surface_crossing_screen import screen
p=np.array([[0,0,0],[1,0,0],[0,1,0],[.1,.1,.2],[.7,.1,.2],[.1,.7,.2]],float)
faces=np.array([[0,1,2],[3,4,5]])
s=SimpleNamespace(free=np.arange(6))
n,gap=reference_axis(p[:3],p[3:])
assert abs(gap-.2)<1e-12 and np.allclose(n,[0,0,1])
R=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
q=p@R.T+[1,2,3]
m,_=transported_axis(p[:3],p[3:],q[:3],q[3:])
assert np.allclose(m,n@R.T)
assert not screen(p,faces)[0]
x=p.copy();x[3,2]=-.1
assert screen(x,faces)[0]=={(0,1)}
c=SelfSeparation(s,faces,p,x,{(0,1)},10.)
assert c.fixedConflictCount==0 and c.A.shape==(9,18)
assert c.violation(x).max()>0
E,g=c.evaluate(x);assert E>0
rng=np.random.default_rng(76);v=rng.normal(size=x.shape);h=1e-6
fd=(c.evaluate(x+h*v,False)[0]-c.evaluate(x-h*v,False)[0])/(2*h)
assert abs(fd-np.sum(g*v))<1e-8*max(1,abs(fd))
H=c.hessian(x);assert v.ravel()@H@v.ravel()>=-1e-10
assert abs(c.evaluate(x+[3,4,5],False)[0]-E)<1e-12
assert SelfSeparation(SimpleNamespace(free=np.array([],dtype=int)),faces,p,x,{(0,1)},10.).fixedConflictCount>0
print('PASS self separation: reference axis, rigid transport, crossing, gradient, PSD, fixed conflict')
