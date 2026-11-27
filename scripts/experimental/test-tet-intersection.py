import numpy as np
from scipy.spatial.transform import Rotation
from tet_intersection import intersection
A = np.r_[np.zeros((1, 3)), np.eye(3)]
clipped = np.array([[.5, -2, -2], [.5, 4, -2], [.5, -2, 4], [4, -2, -2]])
for factor in [1e-3, 1., 1e3]:
    R = Rotation.from_rotvec([.47, -.19, .31]).as_matrix()
    transform = lambda x: (x @ R.T + [17., -23., 11.]) * factor
    for B, expected in [(A, 1 / 6), (A * .5, 1 / 48), (clipped, 1 / 48), (A + [2, 0, 0], 0), (A + [1, 0, 0], 0)]:
        value, _ = intersection(transform(A), transform(B))
        reverse, _ = intersection(transform(B), transform(A))
        assert abs(value - expected * factor**3) <= 1e-9 * factor**3
        assert abs(value - reverse) <= 1e-9 * factor**3
print('PASS tetrahedral overlap integration: independent contained and half-space-clipped volumes, disjoint and point contact, exchange symmetry, rigid/scale covariance')
left = clipped.copy();left[:, 0] = 1 - left[:, 0]
right_volume, _ = intersection(A, clipped)
left_volume, _ = intersection(A, left)
assert abs(right_volume + left_volume - 1 / 6) < 1e-12
assert abs(left_volume - 7 / 48) < 1e-12
F = np.array([[1.7, .2, -.1], [0., .8, .3], [.1, 0., 1.2]])
value, _ = intersection(A @ F.T, clipped @ F.T)
assert abs(value - np.linalg.det(F) / 48) < 1e-12
print('PASS nonmatching material-cell partition: complementary cuts conserve independent source volume; general positive affine covariance')
from tet_intersection import prepare, intersections, intersection_prepared
rng = np.random.default_rng(73076)
all_tets = [A, clipped, left, A * .5, A + [2, 0, 0]]
for i in range(50):
    matrix = np.eye(3) + rng.normal(size=(3, 3)) * .2
    all_tets.append(A @ matrix.T + rng.normal(size=3) * .6)
prepared = prepare(np.array(all_tets))
values, statuses = intersections(prepared[0], prepared[1:])
pairwise = np.array([intersection_prepared(prepared[0], b)[0] for b in prepared[1:]])
assert np.max(np.abs(values - pairwise)) < 1e-12
assert sum(statuses.values()) == len(values)
print('PASS batched integration: independent analytic cases and randomized general cuts agree with pairwise predicates and volumes')
