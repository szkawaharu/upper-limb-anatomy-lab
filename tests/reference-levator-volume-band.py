"""Independent gradient, multiplier and PSD search-metric checks."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experimental'))
from volumetric_reference import Solid
from local_volume_band import VolumeBand
p=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[.2,.2,.2]],float)
t=np.array([[0,1,2,4],[0,1,4,3],[0,4,2,3],[4,1,2,3]])
s=Solid(p,t,[0,1,2,3])
b=VolumeBand(s,30)
E,g=b.evaluate(p)
assert E==0 and np.max(abs(g))==0 and b.hessian(p).nnz==0
x=p.copy();x[4]+=[.04,-.02,.01]
E,g=b.evaluate(x)
assert E>0 and np.all(np.isfinite(g))
rng=np.random.default_rng(76);v=rng.normal(size=x.shape);v[s.fixed]=0
h=1e-6
difference=(b.evaluate(x+h*v,False)[0]-b.evaluate(x-h*v,False)[0])/(2*h)
assert abs(difference-np.sum(g*v))<1e-7*max(1,abs(difference))
H=b.hessian(x);assert np.max(abs((H-H.T).data),initial=0)<1e-10
assert v.ravel()@H@v.ravel()>=-1e-10
assert np.all(b.update(x)>=0)
assert np.any(b.update(x)>0)
print('PASS local volume band: inactive zero, finite-difference gradient, PSD metric, nonnegative multiplier')
