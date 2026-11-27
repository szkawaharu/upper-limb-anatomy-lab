import numpy as np
from landmark_registration import similarity,bone_similarity,apply
p=np.array([[0,0,0],[1,0,0],[0,1,0],[.1,.2,.3]])
R=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
q=1.4*p@R.T+[3,2,-1];a=similarity(p,q)
assert np.max(np.abs(apply(a,p)-q))<1e-12 and np.linalg.det(a['rotation'])>.999999
# Planar scapular anchors are supported without an arbitrary reflected normal.
b=similarity(p[:3],q[:3]);assert np.max(np.abs(apply(b,p)-q))<1e-12
# A mirrored nonplanar atlas must retain residual instead of being reflected.
mirrored=p.copy();mirrored[:,0]*=-1;c=similarity(p,mirrored)
assert np.linalg.det(c['rotation'])>.999999 and np.max(c['errors'])>.01
try:similarity([[0,0,0],[1,0,0],[2,0,0]],[[0,0,0],[0,1,0],[0,2,0]])
except ValueError:pass
else:raise AssertionError('Collinear anchors accepted')
s0=np.array([0.,0,0]);s1=np.array([0.,-2,0]);lat=np.array([0.,0,1]);t0=apply(a,s0);t1=apply(a,s1)
d=bone_similarity(s0,s1,lat,t0,t1,R@lat);assert np.max(np.abs(apply(d,p)-q))<1e-12
print('PASS atlas registration: proper rigid similarity, planar anchors, reflected-atlas residual, degeneracy rejection, bone length/direction without width rescaling')
