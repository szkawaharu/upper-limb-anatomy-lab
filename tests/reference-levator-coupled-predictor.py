"""Independent synthetic check of current coupled boundary predictor."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/experimental'))
from volumetric_reference import Solid
from volumetric_preconditioner import coupled_boundary_predictor,strain_metric_matrix
p=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[.2,.2,.2]],float)
t=np.array([[0,1,2,4],[0,1,4,3],[0,4,2,3],[4,1,2,3]])
s=Solid(p,t,[0,1,2,3])
target=p[s.fixed]+[.03,-.04,.02]
z,report=coupled_boundary_predictor(s,p,target,1000)
assert np.array_equal(z[s.fixed],target)
assert np.max(abs(z-(p+[.03,-.04,.02])))<1e-10
assert report['linearBackwardError']<1e-10
target=p[s.fixed].copy();target[1,0]+=.01;target[2,1]-=.008
z,report=coupled_boundary_predictor(s,p,target,1000)
H=strain_metric_matrix(s,p,1000,full=True)
v=(z-p).ravel();free=(s.free[:,None]*3+np.arange(3)).ravel()
assert np.max(abs((H@v)[free]))<1e-10
assert np.array_equal(z[s.fixed],target)
assert s.evaluate(z,1000,False)[2].min()>0
print('PASS coupled Dirichlet predictor: rigid translation, fixed values, independent free-block residual')
