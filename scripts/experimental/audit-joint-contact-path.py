"""Independent full saved-surface audit of an actual joint-angle path.
Does not certify anatomical attachment data or global self-collision.
"""
import sys
import json
import hashlib
from pathlib import Path
import numpy as np
from material_pose import boundary_targets

p = Path(sys.argv[1])
run = json.loads((p / 'path.json').read_text())
a = np.load(p / 'cage.npz')
s = np.load(p / 'conforming.npz')
rest, tets, faces = a['nodes'], a['tets'], s['faces']
scale, center = float(a['worldUnitsPerMm']), a['centerWorld']
classes = s['classes']
fixed = np.flatnonzero(classes >= 0)
edges = lambda v: np.stack([v[tets[:, k]] - v[tets[:, 0]] for k in [1, 2, 3]], axis=2)
reference_inverse = np.linalg.inv(edges(rest))
source_cells = np.unique(s['elements'])
source = s['points']
reference_normal = np.cross(source[faces[:, 1]] - source[faces[:, 0]], source[faces[:, 2]] - source[faces[:, 0]])
volume = lambda v: float(np.einsum('ti,ti->t', v[faces[:, 0]], np.cross(v[faces[:, 1]], v[faces[:, 2]])).sum() / 6)
rows = []
files = sorted(f for f in p.glob('*-reference-*.npz') if not f.name.endswith('.tmp.npz'))
assert len(files) == len(run['rows'])
for f in files:
    meta = json.loads(f.with_suffix('.json').read_text())
    assert hashlib.sha256(f.read_bytes()).hexdigest() == meta['snapshotSha256']
    assert meta['progress'] == 1 and meta['allStagesConverged'] and meta['contactReached']
    data = np.load(f)
    nodes = data['nodes']
    assert np.isfinite(nodes).all()
    F = edges(nodes) @ reference_inverse
    J = np.linalg.det(F)
    assert np.min(J) > 0
    points = np.einsum('ni,nij->nj', s['embeddingBarycentric'], nodes[s['embeddingVertices']])
    actual_normal = np.cross(points[faces[:, 1]] - points[faces[:, 0]], points[faces[:, 2]] - points[faces[:, 0]])
    selected = F[s['elements']]
    expected = np.linalg.det(selected)[:, None] * np.einsum('tij,tj->ti', np.linalg.inv(selected).transpose(0, 2, 1), reference_normal)
    opposed = np.einsum('ti,ti->t', actual_normal, expected) < 0
    pins = boundary_targets(rest, fixed, classes, meta['target'], center, scale)
    # GPU buffer quantization of control nodes, independent of original surface
    # point order. The runtime GPU deformation solver has not been ported.
    quantized = ((nodes * scale + center).astype(np.float32).astype(float) - center) / scale
    float32_J = np.linalg.det(edges(quantized) @ reference_inverse)
    rows.append({'file': f.name, 'snapshotSha256': meta['snapshotSha256'], 'elev': meta['target']['elev'],
                 'sourceVolumeRatio': volume(points) / volume(source),
                 'attachmentMaxMm': float(np.linalg.norm(nodes[fixed] - pins, axis=1).max(initial=0)),
                 'cageJRange': [float(J.min()), float(J.max())],
                 'surfaceAdjacentJRange': [float(J[source_cells].min()), float(J[source_cells].max())],
                 'surfaceAdjacentCellCount': len(source_cells),
                 'sourceOpposedFaces': int(opposed.sum()), 'cageNonPositive': int((J <= 0).sum()),
                 'float32StorageCageNonPositive': int((float32_J <= 0).sum()),
                 'float32StorageJMin': float(float32_J.min()),
                 'maxDepthIncreaseMm': max(h['maxDepthIncreaseMm'] for h in meta['history'][-1:]),
                 'actualDistanceResidualMm': meta['actualDistanceResidualMm']})
result = {'experimental': True, 'anatomicalAcceptance': False, 'productionEngineAcceptance': False,
          'requestedPathCompleted': run['completed'], 'currentAngle': run['currentAngle'],
          'rows': rows, 'limitations': ['Surface-adjacent local J outliers reported, not hidden in background cell averages',
          'Bone depth copied from final actual-distance solve; independent world-space bone audit is separate',
          'Float32 storage is not a new GPU solver consistency test', 'No global self-collision or mobile FPS certification']}
(p / 'independent-audit.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
