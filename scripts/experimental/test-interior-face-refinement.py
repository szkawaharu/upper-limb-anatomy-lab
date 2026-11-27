import numpy as np
from volumetric_reference import Solid
from interior_face_refinement import refine, metrics, select_faces
from cage_solver import minimize

# Two tetrahedra sharing a face in z=0, positive original orientation.
rest=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,-1]])
tets=np.array([[0,1,2,3],[0,2,1,4]])
solid=Solid(rest,tets,np.array([0,3,4]))
x=rest@np.array([[1.1,.1,0],[0,.9,.1],[0,0,1.2]]).T
faces=select_faces(solid,x,1)
assert faces==[(0,1,2)]
fine,y,parent,info=refine(solid,x,faces)
assert info['exteriorFacesUnchanged'] and info['incidenceFluxRank']==1
assert len(fine.rest)==8 and len(fine.tets)==12
m=metrics(solid,fine,y,parent)
expected=np.linalg.det(solid.edges(x)@solid.inverse)
assert np.allclose(m['parentJ'],expected,atol=1e-12)
assert np.allclose(fine.evaluate(y,30,False)[0],solid.evaluate(x,30,False)[0],atol=1e-12)
# The shared face centre supplies volume exchange even while the original
# boundary positions remain fixed; parent centre alone does not.
z=y.copy();z[5,2]+=.02
after=metrics(solid,fine,z,parent)
difference=(after['parentJ']-m['parentJ'])*solid.volumes
assert abs(difference.sum())<1e-12 and abs(difference[0])>1e-4
assert np.array_equal(z[:len(x)],x)
z=y.copy();z[6,2]+=.02
assert np.allclose(metrics(solid,fine,z,parent)['parentJ'],expected,atol=1e-12)
solution,report=minimize(fine,y,bulk=30,iterations=80)
assert report['converged'],report
assert np.array_equal(solution[solid.fixed],x[solid.fixed])
assert metrics(solid,fine,solution,parent)['allChildJRange'][0]>0
for f in [(0,1,3),(0,2,4)]:
    try:refine(solid,x,[f])
    except ValueError:pass
    else:raise AssertionError('Exterior face accepted')
bad=Solid(rest,tets,np.array([0,3,4]),np.array([.6,.4]))
try:refine(bad,x,faces)
except ValueError:pass
else:raise AssertionError('Nonuniform integration silently inherited')
print('PASS interior face refinement: exact affine transfer, boundary topology, independent parent volume flux, centre-only counterexample, positive-J equilibrium and cut-cell refusal')
