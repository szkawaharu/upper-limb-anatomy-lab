"""Shared actual-bone contact for an embedded source material surface.
Queries are exact mesh occupancy plus point/triangle distance, not an envelope.
Only sampled material points are constrained; global self-collision is separate.
"""
from pathlib import Path
import hashlib
import numpy as np
from contact_potential import PlaneContact, AugmentedSolid
from material_pose import material_to_bone


class FixedContactConflict(ValueError):
    def __init__(self, detail):
        super().__init__('Declared fixed boundary contradicts actual bone contact')
        self.detail = detail


def mesh_distances(mesh, points):
    import trimesh
    inside = mesh.contains(points)
    depth = np.zeros(len(points))
    normal = np.zeros_like(points)
    indices = np.flatnonzero(inside)
    for start in range(0, len(indices), 256):
        use = indices[start:start + 256]
        near, d, face = trimesh.proximity.closest_point(mesh, points[use])
        depth[use] = d
        normal[use] = (near - points[use]) / np.maximum(d[:, None], 1e-30)
        zero = d < 1e-12
        normal[use[zero]] = mesh.face_normals[face[zero]]
    return depth, normal, inside


def atlas_bones(coverage, calibration, center, scale, cache='/tmp/upper-limb-atlas'):
    import trimesh
    assert scale == calibration['scale']
    axes = np.array([[-1, 0, 0], [0, 0, 1], [0, 1, 0]])
    offset = np.asarray(calibration['landmarks']['SC']) / scale - axes @ np.array([-14.91024, -129.66467, 1341.37453])
    result = []
    for entry in coverage['bones']:
        path = Path(cache) / entry['file']
        mesh = trimesh.load_mesh(path, process=True)
        mesh.vertices = np.asarray(mesh.vertices) @ axes.T + offset - center / scale
        if not mesh.is_watertight:
            raise ValueError('Contact bone shell is not closed: ' + entry['file'])
        result.append({**entry, 'query': lambda points, mesh=mesh: mesh_distances(mesh, points),
                       'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    return result


class ActualBoneContact:
    def __init__(self, surface, bones, center, scale):
        self.source = np.asarray(surface['points'])
        self.vertices = np.asarray(surface['embeddingVertices'])
        self.bary = np.asarray(surface['embeddingBarycentric'])
        self.classes = np.asarray(surface['classes'])
        self.bones = bones
        self.center, self.scale = np.asarray(center), scale
        self.baseline = [bone['query'](self.source)[0] for bone in bones]
        faces = np.asarray(surface['faces'])
        tri = self.source[faces]
        area = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2
        self.area_weights = np.zeros(len(self.source))
        for k in range(3):
            np.add.at(self.area_weights, faces[:, k], area / 3)
        self.area_weights /= self.area_weights.sum()
        edges = np.unique(np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1), axis=0)
        spacing = np.linalg.norm(self.source[edges[:, 1]] - self.source[edges[:, 0]], axis=1)
        self.initial_penalty = 1 / float(np.median(spacing))**2

    def surface(self, nodes):
        return np.einsum('ni,nij->nj', self.bary, nodes[self.vertices])

    def scan(self, nodes, target, make_planes=False):
        points = self.surface(nodes)
        rows, planes = [], []
        for bone, baseline in zip(self.bones, self.baseline):
            R, t = material_to_bone(target, bone['key'], self.center, self.scale)
            depth, normal, inside = bone['query'](points @ R + t)
            increase = depth - baseline
            worst = int(increase.argmax())
            rows.append({'bone': bone['file'], 'maxBaselineDepthMm': float(baseline.max()),
                         'maxDepthMm': float(depth.max()), 'maxDepthIncreaseMm': float(increase.max()),
                         'worstSurfaceVertex': worst, 'insidePoints': int(inside.sum()),
                         'newInsidePoints': int((inside & (baseline == 0)).sum())})
            if not make_planes:
                continue
            use = np.flatnonzero(inside)
            n = normal[use] @ R.T
            offsets = np.einsum('ni,ni->n', n, points[use]) + increase[use]
            immovable = np.all((self.classes[self.vertices[use]] >= 0) | (self.bary[use] <= 1e-9), axis=1)
            conflict = immovable & (increase[use] > 1e-5)
            if conflict.any():
                bad = use[conflict]
                raise FixedContactConflict({'bone': bone['file'], 'vertices': bad.tolist(),
                                            'maxDepthIncreaseMm': float(increase[bad].max())})
            planes.append((use, n, offsets))
        return rows, planes

    def augmented(self, solid, nodes, target, penalty):
        rows, planes = self.scan(nodes, target, True)
        use = np.concatenate([p[0] for p in planes])
        normals = np.concatenate([p[1] for p in planes])
        offsets = np.concatenate([p[2] for p in planes])
        potential = PlaneContact(len(nodes), self.vertices[use], self.bary[use], normals, offsets, self.area_weights[use], penalty)
        return AugmentedSolid(solid, [potential]), rows

    def equilibrate(self, solid, initial, target, minimize, bulk=30., tolerance_mm=.1, iterations=150):
        """Use actual mesh distance in every energy/gradient/line-search call.
        The contact Hessian is its positive Gauss-Newton approximation. This
        changes the direction only; true energy, residual and positive-J guards
        decide whether a trial is accepted. Penalty levels never change pins.
        """
        x = initial.copy()
        fixed = x[solid.fixed].copy()
        penalty = self.initial_penalty
        history = []
        for level in range(8):
            potential = ActualDistancePotential(self, len(x), target, penalty)
            augmented = AugmentedSolid(solid, [potential])
            trial, report = minimize(augmented, x, bulk=bulk, iterations=iterations)
            assert np.array_equal(trial[solid.fixed], fixed)
            _, gradient, J = augmented.evaluate(trial, bulk)
            residual = float(np.linalg.norm(solid.factor.solve(gradient[solid.free]), axis=1).max(initial=0))
            rows = potential.state(trial)[1]
            maximum = max(r['maxDepthIncreaseMm'] for r in rows)
            history.append({'penaltyLevel': level, 'penaltyPerMm2': penalty,
                            'actualDistanceResidualMm': residual, 'maxDepthIncreaseMm': maximum,
                            'bones': rows, 'minimization': report})
            if not report['converged'] or residual >= .001:
                return initial.copy(), {'converged': False, 'reason': 'actual-distance-minimization-failed', 'history': history}
            x = trial
            if maximum <= tolerance_mm:
                return x, {'converged': True, 'reason': 'actual-distance-equilibrium-and-contact', 'history': history,
                           'contactReached': True, 'contactToleranceMm': tolerance_mm,
                           'actualDistanceResidualMm': residual, 'cageJQuantiles': np.quantile(J, [0, .01, .5, .99, 1]).tolist()}
            penalty *= 10
        return initial.copy(), {'converged': False, 'reason': 'contact-penalty-limit', 'history': history}


class ActualDistancePotential:
    def __init__(self, contact, node_count, target, penalty):
        self.contact, self.node_count, self.target, self.penalty = contact, node_count, target, penalty
        self.cached_nodes = None

    def state(self, nodes):
        if self.cached_nodes is None or not np.array_equal(nodes, self.cached_nodes):
            rows, planes = self.contact.scan(nodes, self.target, True)
            use = np.concatenate([p[0] for p in planes])
            normals = np.concatenate([p[1] for p in planes])
            offsets = np.concatenate([p[2] for p in planes])
            potential = PlaneContact(self.node_count, self.contact.vertices[use], self.contact.bary[use],
                                     normals, offsets, self.contact.area_weights[use], self.penalty)
            self.cached_nodes = nodes.copy()
            self.cached_value = potential, rows
        return self.cached_value

    def evaluate(self, nodes, gradient=True):
        return self.state(nodes)[0].evaluate(nodes, gradient)

    def hessian(self, nodes):
        # PSD Gauss-Newton metric; omitted signed-distance curvature is not
        # presented as an exact material/contact Hessian.
        return self.state(nodes)[0].hessian(nodes)
