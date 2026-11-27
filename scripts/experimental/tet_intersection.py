"""Convex tetrahedron intersection using edge/face candidates and hull volume.
Offline floating-point integration, not exact rational geometry. Per-tet local
normalization and cached inverses avoid repeated factorizations for shared cells.
"""
import numpy as np
from scipy.spatial import ConvexHull

EDGES = np.array([[i, j] for i in range(4) for j in range(i + 1, 4)])


def prepare(tets):
    points = np.asarray(tets, float)
    scale = np.max(np.ptp(points, axis=1), axis=1)
    if not np.isfinite(points).all() or np.any(scale <= 0):
        raise ValueError('Nonfinite or collapsed tetrahedron')
    D = np.stack([points[:, k] - points[:, 0] for k in [1, 2, 3]], axis=2) / scale[:, None, None]
    inverse = np.linalg.inv(D)
    condition = np.linalg.cond(D)
    volume = np.abs(np.linalg.det(D)) * scale**3 / 6
    return [(p, s, inv, cond, v) for p, s, inv, cond, v in zip(points, scale, inverse, condition, volume)]


def intersection_prepared(first, second):
    A, B = first[0], second[0]
    tolerance = 64 * np.finfo(float).eps * max(first[3], second[3])
    def bary(points, prepared):
        tet, scale, inverse = prepared[:3]
        three = ((points - tet[0]) / scale) @ inverse.T
        return np.c_[1 - three.sum(1), three]
    candidates = []
    for primary, own, other in [(A, first, second), (B, second, first)]:
        coordinates = bary(primary, other)
        if np.any(np.max(coordinates, axis=0) < -tolerance):
            return 0., {'status': 'separated-by-tetrahedral-face', 'barycentricRoundoffScale': tolerance}
        inside = np.all(coordinates >= -tolerance, axis=1)
        if np.all(inside):
            return float(own[4]), {'status': 'contained-tetrahedron', 'barycentricRoundoffScale': tolerance}
        candidates.extend(primary[inside])
        i, j = EDGES.T
        denominator = coordinates[i] - coordinates[j]
        with np.errstate(divide='ignore', invalid='ignore'):
            t = coordinates[i] / denominator
        valid = np.isfinite(t) & (t >= 0) & (t <= 1)
        edge, face = np.nonzero(valid)
        point = primary[i[edge]] + t[edge, face, None] * (primary[j[edge]] - primary[i[edge]])
        candidates.extend(point[np.all(bary(point, other) >= -tolerance, axis=1)])
    points = np.asarray(candidates)
    if len(points) < 4:
        return 0., {'status': 'empty-or-lower-dimensional', 'barycentricRoundoffScale': tolerance}
    scale = np.max(np.ptp(points, axis=0))
    if scale == 0:
        return 0., {'status': 'rank-unresolved-or-contact', 'barycentricRoundoffScale': tolerance}
    normalized = (points - points[0]) / scale
    centered = normalized - normalized.mean(0)
    _, singular, axes = np.linalg.svd(centered, full_matrices=False)
    rank_error = 64 * np.finfo(float).eps * max(centered.shape) * max(singular[0], 1.)
    if singular[-1] <= rank_error:
        return 0., {'status': 'rank-unresolved-or-contact', 'rankRoundoffScale': rank_error, 'barycentricRoundoffScale': tolerance}
    # A rigid PCA chart plus anisotropic normalization prevents Qhull from
    # discarding a small but resolvable third dimension. Undo its determinant
    # in the returned physical volume; this does not thicken the geometry.
    whitened = (centered @ axes.T) / singular
    hull = ConvexHull(whitened)
    return float(hull.volume) * float(singular.prod()) * scale**3, {'status': 'positive-volume', 'barycentricRoundoffScale': tolerance, 'candidatePoints': len(points)}


def intersection(a, b):
    first, second = prepare(np.array([a, b]))
    return intersection_prepared(first, second)


def intersections(first, others):
    """Vectorized separating-face/containment tests, same pairwise predicates.
    Only unresolved edge/face cuts need the convex-hull construction.
    """
    if not others:
        return np.empty(0), {}
    tet, scale, inverse, condition, volume = first
    points = np.stack([o[0] for o in others])
    scales = np.array([o[1] for o in others])
    inverses = np.stack([o[2] for o in others])
    tolerance = 64 * np.finfo(float).eps * np.maximum(condition, [o[3] for o in others])
    three = np.einsum('nvi,nji->nvj', (tet[None] - points[:, 0, None]) / scales[:, None, None], inverses)
    A_in_B = np.concatenate([1 - three.sum(2, keepdims=True), three], axis=2)
    three = ((points - tet[0]) / scale) @ inverse.T
    B_in_A = np.concatenate([1 - three.sum(2, keepdims=True), three], axis=2)
    separated = np.any(A_in_B.max(1) < -tolerance[:, None], axis=1) | np.any(B_in_A.max(1) < -tolerance[:, None], axis=1)
    first_inside = (~separated) & np.all(A_in_B >= -tolerance[:, None, None], axis=(1, 2))
    second_inside = (~separated) & (~first_inside) & np.all(B_in_A >= -tolerance[:, None, None], axis=(1, 2))
    result = np.zeros(len(others))
    result[first_inside] = volume
    result[second_inside] = np.array([o[4] for o in others])[second_inside]
    statuses = {'separated-by-tetrahedral-face': int(separated.sum()), 'contained-tetrahedron': int(first_inside.sum() + second_inside.sum())}
    for i in np.flatnonzero(~(separated | first_inside | second_inside)):
        result[i], report = intersection_prepared(first, others[i])
        status = report['status'];statuses[status] = statuses.get(status, 0) + 1
    return result, statuses
