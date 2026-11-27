import numpy as np
from volumetric_reference import Solid
from local_volume_constraints import LocalPressure, log_volume_jacobian, physical_cells, audit, equilibrate
from actual_bone_contact import ActualBoneContact

rest = np.array([[3.,0,0],[3.,0,1],[3.,1,0],[2.1,.2,.2]])
tets = np.array([[0,1,2,3]])
solid = Solid(rest, tets, np.arange(3))
cells, w = physical_cells(solid, np.array([.15]), np.array([.02]))
assert cells.tolist() == [0]
p = LocalPressure(solid, cells, w, np.array([.7]), 30.)
x = rest.copy(); x[3,0] += .2
c,A = log_volume_jacobian(solid, x, cells)
for k in range(x.size):
    plus,minus=x.copy(),x.copy()
    plus.flat[k]+=1e-6;minus.flat[k]-=1e-6
    numeric=(log_volume_jacobian(solid,plus,cells)[0]-log_volume_jacobian(solid,minus,cells)[0])/2e-6
    assert np.max(np.abs(numeric-A.toarray()[:,k])) < 1e-8
    numeric=(p.evaluate(plus,False)[0]-p.evaluate(minus,False)[0])/2e-6
    assert abs(numeric-p.evaluate(x)[1].flat[k]) < 1e-7
assert np.linalg.eigvalsh(p.hessian(x).toarray()).min() > -1e-10
assert audit(solid,x,cells)['linearizedResidualMax'] < 1e-10
try:
    physical_cells(solid,np.array([0.]),np.array([1.]))
except ValueError:
    pass
else:
    raise AssertionError('Invalid material integration accepted')
# Independent analytic sphere; actual-distance contact is evaluated at each trial.
def sphere(points):
    d=np.linalg.norm(points,axis=1)
    return np.maximum(2-d,0), points/np.maximum(d[:,None],1e-30), d<2
surface={'points':rest,'faces':np.array([[1,2,3],[0,3,2],[0,1,3],[0,2,1]]),
         'embeddingVertices':np.tile(np.arange(4),(4,1)),
         'embeddingBarycentric':np.eye(4),'classes':np.array([0,0,0,-1])}
contact=ActualBoneContact(surface,[{'key':'b','file':'sphere','query':sphere}],np.zeros(3),1.)
identity={'q':[0,0,0,1],'p':[0,0,0],'rest':[0,0,0]}
target={'a':identity,'b':{**identity,'p':[.3,0,0]}}
sol,report=equilibrate(solid,x,contact,target,np.array([.15]),np.array([.02]),
                      max_log_error=.01,outer_steps=30)
assert report['converged'],report
assert not report['accepted']
assert np.array_equal(sol[:3],x[:3])
assert abs(np.log(solid.evaluate(sol)[2][0]))<=.01
assert contact.scan(sol,target)[0][0]['maxDepthIncreaseMm']<=.1
# An intentionally impossible outer iteration budget must not leak a trial.
failed,report=equilibrate(solid,x,contact,target,np.array([.15]),np.array([.02]),
                         max_log_error=1e-14,outer_steps=1)
assert not report['converged'] and np.array_equal(failed,x)
print('PASS local volume derivatives, PSD metric, linear feasibility, actual-distance coupled solve and rollback')

# A tiny physical intersection must not disappear into the background, and
# an all-fixed deformed element cannot be silently freed or satisfied by
# pressure in a different element.
from types import SimpleNamespace
two = np.concatenate([rest, rest+[10.,0,0]], axis=0)
volumes = np.array([1e-14, 1.])
background = np.array([.02,.02])
weights = volumes/volumes.sum()*(1-background.sum())+background
solid2 = Solid(two, np.array([[0,1,2,3],[4,5,6,7]]), np.arange(7), weights)
cells2, w2 = physical_cells(solid2, volumes, background)
assert cells2.tolist() == [0,1] and 0 < w2[0] < 1e-13
deformed = two.copy(); deformed[3,0]+=.1
a = audit(solid2,deformed,cells2)
assert a['immovableCellIds']==[0] and a['constraintRankDeficiencyLowerBound']==1
assert a['immovableMaxAbsLogJ']>.01
failed, report = equilibrate(solid2,deformed,SimpleNamespace(initial_penalty=1.),target,volumes,background,max_log_error=.01)
assert report['reason']=='fixed-local-volume-conflict' and np.array_equal(failed,deformed)
print('PASS tiny physical cell retained and fully fixed incompatible volume rejected')
