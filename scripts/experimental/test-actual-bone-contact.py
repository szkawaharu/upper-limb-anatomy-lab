import numpy as np
from actual_bone_contact import ActualBoneContact, FixedContactConflict
from material_pose import rigid
from scipy.spatial.transform import Rotation

# Analytic closed sphere: independent test oracle, no mesh library dependency.
def sphere(points):
    distance = np.linalg.norm(points, axis=1)
    inside = distance < 2
    depth = np.maximum(2 - distance, 0)
    normals = points / np.maximum(distance[:, None], 1e-30)
    return depth, normals, inside
nodes = np.array([[2.1, 0., 0.], [2.1, .1, 0], [2.1, 0, .1], [2.2, 0, 0]])
ids = np.tile(np.arange(4), (3, 1))
bary = np.eye(4)[:3]
surface = {'points': nodes[:3], 'faces': np.array([[0, 1, 2]]), 'embeddingVertices': ids, 'embeddingBarycentric': bary, 'classes': np.array([-1] * 4)}
bones = [{'key': 'b', 'file': 'analytic-sphere', 'query': sphere}]
contact = ActualBoneContact(surface, bones, np.zeros(3), 1.)
identity = {'q': [0, 0, 0, 1], 'p': [0, 0, 0], 'rest': [0, 0, 0]}
target = {'a': identity, 'b': identity}
x = nodes.copy();x[:, 0] -= .3
rows, planes = contact.scan(x, target, True)
assert abs(rows[0]['maxDepthIncreaseMm'] - .2) < 1e-12
use, normal, offsets = planes[0]
assert np.max(np.abs(offsets - np.einsum('ni,ni->n', normal, x[use]) - (2 - np.linalg.norm(x[use], axis=1)))) < 1e-12
# Independent rigidly transformed geometry: material bone query is identical.
R = Rotation.from_rotvec([.31, -.27, .2]);shift = np.array([4., 2., -1.])
frame = {'q': R.as_quat(), 'p': shift, 'rest': [0, 0, 0]}
rows2, _ = contact.scan(x, {'a': frame, 'b': frame}, False)
assert abs(rows2[0]['maxDepthIncreaseMm'] - rows[0]['maxDepthIncreaseMm']) < 1e-12
surface['classes'] = np.array([0, 0, 0, -1])
fixed = ActualBoneContact(surface, bones, np.zeros(3), 1.)
try:
    fixed.scan(x, target, True)
except FixedContactConflict as error:
    assert error.detail['vertices'] == [0, 1, 2]
else:
    raise AssertionError('Fixed region was silently displaced or unpinned')
# Baseline overlap is retained pointwise, not converted into a muscle-wide offset.
surface['points'] = x[:3].copy();surface['classes'] = np.array([-1] * 4)
baseline = ActualBoneContact(surface, bones, np.zeros(3), 1.)
assert abs(baseline.scan(x, target)[0][0]['maxDepthIncreaseMm']) < 1e-12
print('PASS actual-bone contact: analytic sphere depth and plane oracle, rigid-frame covariance, pointwise baseline allowance, contradictory fixed region rejected without unpinning')

# The refreshed tangent potential has the exact first derivative of sampled
# squared actual sphere depth at the current geometry, independently differenced.
from scipy.sparse import eye
class InertMaterial:
    fixed = np.array([], dtype=int)
    free = np.arange(4)
    factor = eye(4, format='csc')
    def evaluate(self, x, bulk=30., gradient=True):
        return 0., np.zeros_like(x) if gradient else None, np.ones(1)
    class Identity:
        def solve(self, rhs):
            return rhs
    factor = Identity()
base = InertMaterial()
augmented, _ = contact.augmented(base, x, target, 100.)
_, g, _ = augmented.evaluate(x)
def exact_energy(y):
    depth = sphere(contact.surface(y))[0]
    return 50. * np.dot(contact.area_weights, depth**2)
for k in range(x.size):
    h = 1e-6
    plus, minus = x.copy(), x.copy()
    plus.flat[k] += h;minus.flat[k] -= h
    derivative = (exact_energy(plus) - exact_energy(minus)) / (2 * h)
    assert abs(derivative - g.flat[k]) < 1e-6
# A stale local minimizer's 'converged' flag cannot certify actual contact or
# equilibrium. Leave geometry untouched to simulate this invalid solver claim.
def falsely_converged(solid, initial, **kwargs):
    return initial.copy(), {'converged': True}
rejected, report = contact.equilibrate(base, x, target, falsely_converged)
assert not report['converged'] and report['reason'] == 'actual-distance-minimization-failed'
assert np.array_equal(rejected, x)
print('PASS refreshed actual-distance stopping: independent nonlinear energy derivative; stale tangent convergence rejected with rollback')

from actual_bone_contact import ActualDistancePotential
actual = ActualDistancePotential(contact, len(x), target, 100.)
for k in range(x.size):
    h = 1e-6;plus, minus = x.copy(), x.copy();plus.flat[k] += h;minus.flat[k] -= h
    derivative = (actual.evaluate(plus, False)[0] - actual.evaluate(minus, False)[0]) / (2 * h)
    assert abs(derivative - actual.evaluate(x)[1].flat[k]) < 1e-6
    assert abs(actual.evaluate(plus, False)[0] - exact_energy(plus)) < 1e-12
print('PASS actual-distance objective: independently differenced energy/gradient at each trial; cached planes never freeze the line-search geometry')

# Full material/contact solve on an independently defined tetrahedron. Three
# fixed nodes remain outside a translated sphere; its free tip must slide out.
from volumetric_reference import Solid
from cage_solver import minimize
rest = np.array([[3., 0, 0], [3., 0, 1], [3., 1, 0], [2.1, .2, .2]])
assert np.linalg.det(np.stack([rest[i] - rest[0] for i in [1, 2, 3]], axis=1)) > 0
solid = Solid(rest, np.array([[0, 1, 2, 3]]), np.arange(3))
source = {'points': rest, 'faces': np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]]),
          'embeddingVertices': np.tile(np.arange(4), (4, 1)), 'embeddingBarycentric': np.eye(4),
          'classes': np.array([0, 0, 0, -1])}
coupled = ActualBoneContact(source, bones, np.zeros(3), 1.)
shifted = {'a': identity, 'b': {**identity, 'p': [.3, 0, 0]}}
solution, result = coupled.equilibrate(solid, rest, shifted, minimize)
assert result['converged'] and result['contactReached'], result
assert np.array_equal(solution[:3], rest[:3])
assert np.isfinite(solution).all() and np.min(solid.evaluate(solution)[2]) > 0
assert coupled.scan(solution, shifted)[0][0]['maxDepthIncreaseMm'] <= .1
assert result['actualDistanceResidualMm'] < .001
print('PASS coupled curved contact/material equilibrium: exact pins, positive J, independently queried sphere clearance and true residual')
