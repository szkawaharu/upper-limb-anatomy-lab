"""Independent finite-difference and rigid-motion checks for StretchBand."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experimental'))
from volumetric_reference import Solid
from local_stretch_band import StretchBand
p=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[.2,.2,.2]],float)
t=np.array([[0,1,2,4],[0,1,4,3],[0,4,2,3],[4,1,2,3]])
s=Solid(p,t,[0,1,2,3])
b=StretchBand(s,10)
E,g=b.evaluate(p);assert E<1e-20 and np.max(abs(g))<1e-10
R=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
E,g=b.evaluate(p@R.T+[2,3,4]);assert E<1e-20 and np.max(abs(g))<1e-10
x=p@np.diag([3.,1.,1/3.]);x[4]+=[.01,-.008,.002]
assert s.evaluate(x,1000,False)[2].min()>0
E,g=b.evaluate(x);assert E>0 and np.all(np.isfinite(g))
rng=np.random.default_rng(76);v=rng.normal(size=x.shape);h=1e-6
fd=(b.evaluate(x+h*v,False)[0]-b.evaluate(x-h*v,False)[0])/(2*h)
assert abs(fd-np.sum(g*v))<1e-7*max(1,abs(fd))
H=b.hessian(x);assert np.max(abs((H-H.T).data),initial=0)<1e-9
assert v.ravel()@H@v.ravel()>=-1e-9
assert abs(b.evaluate(x@R.T+[2,3,4],False)[0]-E)<1e-10
print('PASS stretch band: rigid invariance, gradient, positive search metric')
