import numpy as np
from joint_path import follow
from material_pose import material_to_bone, material_to_world, boundary_targets, rigid
from scipy.spatial.transform import Rotation

rng = np.random.default_rng(73076)
center = rng.normal(size=3)
scale = .00555
frame = lambda: {'q': Rotation.random(random_state=rng).as_quat(), 'p': rng.normal(size=3), 'rest': rng.normal(size=3)}
target = {'a': frame(), 'b': frame()}
points = rng.normal(size=(13, 3)) * 50
RA, ta = rigid(target['a'])
RB, tb = rigid(target['b'])
world = material_to_world(points, target, center, scale)
R, t = material_to_bone(target, 'b', center, scale)
query = points @ R + t
independent = ((world - tb) @ RB - center) / scale
assert np.max(np.abs(query - independent)) < 1e-10
classes = np.array([0, 1] * 6 + [-1])
fixed = np.flatnonzero(classes >= 0)
boundary = boundary_targets(points, fixed, classes, target, center, scale)
posed = material_to_world(boundary, target, center, scale)
expected = (points[fixed] * scale + center) @ RB.T + tb
assert np.max(np.abs(posed[classes[fixed] == 1] - expected[classes[fixed] == 1])) < 1e-12
assert np.array_equal(boundary[classes[fixed] == 0], points[fixed][classes[fixed] == 0])

observed = []
def solve(x, previous, next_index):
    assert x[0] == previous, 'Rejected or stale trial leaked into the accepted state'
    x[0] = next_index
    return x, {'converged': abs(next_index - previous) <= 2}
x, report = follow(np.array([0.]), 0, [8, 0], solve, lambda x, i, a, r: observed.append(i))
assert report['completed'] and x[0] == 0 and 8 in observed and observed[-1] == 0
assert all(0 < abs(r['index'] - r['fromIndex']) <= 2 for r in report['accepted'])
assert report['rejected']
def blocked(x, previous, target):
    x[:] = 999
    return x, {'converged': False}
x, report = follow(np.array([7.]), 7, [10], blocked)
assert not report['completed'] and report['index'] == 7 and x[0] == 7
print('PASS actual joint path: independent world/material/bone transforms, exact broad boundary, adaptive actual-grid subdivision, reverse traversal, rejected-state rollback and explicit incomplete endpoint')

# Multiple cervical/trunk/scapular attachment regions use their own final bones.
# Synthetic frames check the contract, not anatomical cervical kinematics.
keys = {0: 'trunk', 1: 'C0', 2: 'C4', 3: 'C7', 4: 'scapula'}
multi = {'a': frame(), **{key: frame() for key in keys.values()}}
labels = np.arange(len(points)) % len(keys);fixed_multi = np.arange(len(points))
bound = boundary_targets(points, fixed_multi, labels, multi, center, scale, keys)
world = material_to_world(bound, multi, center, scale)
for label, key in keys.items():
    R, t = rigid(multi[key]);use = labels == label
    expected = (points[use] * scale + center) @ R.T + t
    assert np.max(np.abs(world[use] - expected)) < 1e-12
mirror = np.diag([-1., 1, 1])
left = {}
for key, f in multi.items():
    R, _ = rigid(f)
    left[key] = {'q': Rotation.from_matrix(mirror @ R @ mirror).as_quat(),
                 'p': mirror @ f['p'], 'rest': mirror @ f['rest']}
mirrored = boundary_targets(points @ mirror, fixed_multi, labels, left, mirror @ center, scale, keys)
assert np.max(np.abs(mirrored - bound @ mirror)) < 1e-10
try:
    boundary_targets(points, fixed_multi, labels, multi, center, scale)
except ValueError:
    pass
else:
    raise AssertionError('An unknown multi-bone attachment silently followed bone a')
print('PASS multiple-bone attachment transport: independent trunk/C0/C4/C7/scapula frames, mirrored side covariance, unregistered classes rejected; synthetic, not anatomical cervical validation')
