"""Offline actual-joint continuation with bone contact at every accepted pose.
Uses frozen final rig frames on a dense angle grid. Endpoint bone SLERP is not
used to approximate intermediate joint angles. This is not a real-time backend.
"""
import sys
import json
import hashlib
import shutil
import time
from pathlib import Path
import numpy as np
from actual_bone_contact import ActualBoneContact, FixedContactConflict, atlas_bones
from bone_sets import bone_set, model_id
from material_pose import boundary_targets, material_to_world
from joint_path import follow
from volumetric_reference import Solid
from cage_solver import minimize
from snapshot_names import snapshot_name


root = Path(sys.argv[1])
initial_path = Path(sys.argv[2])
out = Path(sys.argv[3])
requested = json.loads(sys.argv[4])
plane = float(sys.argv[5]) if len(sys.argv) > 5 else 0.
clock = time.time()
a = np.load(root / 'cage.npz')
s = np.load(root / 'conforming.npz')
initial = np.load(initial_path)
initial_report = json.loads(initial_path.with_suffix('.json').read_text())
assert initial_report['progress'] == 1 and initial_report['allStagesConverged']
assert initial_report['contactReached']
assert hashlib.sha256(initial_path.read_bytes()).hexdigest() == initial_report['snapshotSha256']
assert initial_report['target']['plane'] == plane
geometry_hashes = dict(initial_report.get('inputSha256', {}))
for path, sha in initial_report.get('inputsSha256', {}).items():
    name = Path(path).name
    if name in ['cage.npz', 'conforming.npz']:
        if name in geometry_hashes:
            assert geometry_hashes[name] == sha
        geometry_hashes[name] = sha
for name in ['cage.npz', 'conforming.npz']:
    assert hashlib.sha256((root / name).read_bytes()).hexdigest() == geometry_hashes[name]
assert requested[0] == initial_report['target']['elev']
assert requested and all(isinstance(v, (float, int)) and np.isfinite(v) and 0 <= v <= 180 for v in requested)
out.mkdir(parents=True, exist_ok=True)
target_file = out / 'targets.json'
targets = json.loads(target_file.read_text())
assert targets['id'] == model_id(initial_path, initial_report)
coverage = bone_set(targets['id'])
poses = sorted([p for p in targets['poses'] if p['plane'] == plane], key=lambda p: p['elev'])
angles = np.array([p['elev'] for p in poses])
assert len(angles) > 1 and np.allclose(np.diff(angles), .375, rtol=0, atol=1e-12)
index_for = lambda value: int(np.flatnonzero(np.abs(angles - value) < 1e-10)[0])
start_index = index_for(requested[0])
for key in ['a', 'b']:
    for value in ['q', 'p', 'rest']:
        assert np.max(np.abs(np.asarray(poses[start_index][key][value]) - initial_report['target'][key][value])) < 1e-9
rest, classes = a['nodes'], s['classes']
fixed = np.flatnonzero(classes >= 0)
solid = Solid(rest, a['tets'], fixed, element_weights=a['weights'])
center, scale = a['centerWorld'], float(a['worldUnitsPerMm'])
bulk = initial_report['bulk']
cal_path = Path('docs/validation/v076/source-bones/browser-numeric/rig-calibration.json')
cal = json.loads(cal_path.read_text())
contact = ActualBoneContact(s, atlas_bones(coverage, cal, center, scale), center, scale)
faces = s['faces']
source = s['points']
volume = lambda y: float(np.einsum('ti,ti->t', y[faces[:, 0]], np.cross(y[faces[:, 1]], y[faces[:, 2]])).sum() / 6)
source_volume = volume(source)
reference_normal = np.cross(source[faces[:, 1]] - source[faces[:, 0]], source[faces[:, 2]] - source[faces[:, 0]])
for name in ['cage.npz', 'conforming.npz', 'preparation.json']:
    shutil.copy2(root / name, out / name)
configuration_hashes = {name: hashlib.sha256((out / name).read_bytes()).hexdigest() for name in ['cage.npz', 'conforming.npz', 'targets.json', 'preparation.json']}
code_files = ['run-joint-contact-path.py', 'joint_path.py', 'material_pose.py', 'actual_bone_contact.py', 'contact_potential.py', 'bone_sets.py', 'contact-bone-sets.json', 'cage_solver.py', 'projected_solid_metric.py', 'constitutive_metric.py', 'volumetric_reference.py', 'snapshot_names.py']
code_hashes = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in code_files}
(out / 'code-snapshot').mkdir(exist_ok=True)
for name in code_files:
    shutil.copy2(Path(__file__).with_name(name), out / 'code-snapshot' / name)
input_hashes = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in [initial_path, initial_path.with_suffix('.json'), target_file, cal_path, root / 'cage.npz', root / 'conforming.npz']}
initial_nodes = initial['nodes'].copy()
first_world = material_to_world(contact.surface(initial_nodes), poses[start_index], center, scale)
seen = {start_index: first_world.copy()}
previous_world = first_world.copy()
rows = []


def solve_at(x, previous, index):
    target = poses[index]
    pins = boundary_targets(rest, fixed, classes, target, center, scale)
    start = solid.predict_boundary(x, pins)
    if not np.isfinite(solid.evaluate(start, bulk, gradient=False)[0]):
        return x, {'converged': False, 'reason': 'nonpositive-J-boundary-predictor'}
    try:
        trial, report = contact.equilibrate(solid, start, target, minimize, bulk=bulk)
    except FixedContactConflict as error:
        return x, {'converged': False, 'reason': str(error), 'fixedConflict': error.detail}
    print(json.dumps({'fromAngle': angles[previous], 'toAngle': angles[index], 'seconds': time.time() - clock,
                      **{k: v for k, v in report.items() if k != 'history'}}), flush=True)
    return trial, report


def save(nodes, index, accepted, rejected):
    global previous_world
    target = poses[index]
    points = contact.surface(nodes)
    world = material_to_world(points, target, center, scale)
    F = (solid.edges(nodes) @ solid.inverse)[s['elements']]
    expected_normal = np.linalg.det(F)[:, None] * np.einsum('tij,tj->ti', np.linalg.inv(F).transpose(0, 2, 1), reference_normal)
    actual_normal = np.cross(points[faces[:, 1]] - points[faces[:, 0]], points[faces[:, 2]] - points[faces[:, 0]])
    prescribed = boundary_targets(rest, fixed, classes, target, center, scale)
    row = {'elev': target['elev'], 'plane': plane, 'sequence': len(rows),
           'surfaceVolumeRatio': volume(points) / source_volume,
           'sourceOpposedFaces': int((np.einsum('ti,ti->t', expected_normal, actual_normal) < 0).sum()),
           'attachmentMaxMm': float(np.linalg.norm(nodes[fixed] - prescribed, axis=1).max(initial=0)),
           'maxStepWorldMm': float(np.linalg.norm(world - previous_world, axis=1).max(initial=0) / scale),
           'returnToPreviousSamePoseMaxMm': None if index not in seen else float(np.linalg.norm(world - seen[index], axis=1).max(initial=0) / scale),
           'seconds': time.time() - clock, **accepted[-1]}
    seen.setdefault(index, world.copy())
    previous_world = world.copy()
    rows.append(row)
    key = str(len(rows)).zfill(3) + '-' + snapshot_name(target['elev'], plane, bulk)
    path = out / key
    temporary = path.with_suffix('.tmp.npz')
    np.savez_compressed(temporary, nodes=nodes, referenceNodes=rest, cageTets=a['tets'], sourceWorld=world,
                        referenceWorld=source * scale + center, sourceFaces=faces, sourceWeights=s['weights'],
                        embeddingVertices=s['embeddingVertices'], embeddingBarycentric=s['embeddingBarycentric'],
                        triangleElements=s['elements'], progress=1.)
    sha = hashlib.sha256(temporary.read_bytes()).hexdigest()
    temporary.replace(path.with_suffix('.npz'))
    report = {'experimental': True, 'accepted': False, 'modelId': targets['id'], 'contactModelId': targets['id'],
              'target': target, 'progress': 1., 'bulk': bulk, 'allStagesConverged': True,
              'contactReached': True, 'boneContactEvaluated': True, 'selfCollisionEvaluated': False,
              'contactBoneCoverage': coverage, 'sourceVertices': len(source), 'originalSourceVertices': len(a['sourcePoints']),
              'cageNodes': len(rest), 'snapshotSha256': sha, 'codeSha256': code_hashes, 'inputsSha256': input_hashes, 'inputSha256': configuration_hashes,
              'activation': None, 'force': None, **row}
    temporary = path.with_suffix('.tmp.json')
    temporary.write_text(json.dumps(report, indent=2) + '\n')
    temporary.replace(path.with_suffix('.json'))
    summary = {'experimental': True, 'accepted': False, 'completed': False, 'currentAngle': angles[index],
               'requestedAngles': requested, 'plane': plane, 'rows': rows, 'rejected': rejected,
               'initialSnapshot': str(initial_path), 'inputsSha256': input_hashes, 'codeSha256': code_hashes,
               'boneCoverage': coverage, 'actualJointGridStepDegrees': .375, 'seconds': time.time() - clock}
    temporary = out / 'path.tmp.json'
    temporary.write_text(json.dumps(summary, indent=2) + '\n')
    temporary.replace(out / 'path.json')


x, result = follow(initial_nodes, start_index, [index_for(v) for v in requested[1:]], solve_at, save)
summary_path = out / 'path.json'
summary = json.loads(summary_path.read_text()) if summary_path.exists() else {'rows': [], 'requestedAngles': requested}
summary.update({'completed': result['completed'], 'currentAngle': angles[result['index']], 'rejected': result['rejected'],
                'seconds': time.time() - clock, 'experimental': True, 'accepted': False})
summary_path.write_text(json.dumps(summary, indent=2) + '\n')
print('FINISHED', json.dumps({k: v for k, v in summary.items() if k not in ['rows', 'rejected', 'inputsSha256', 'codeSha256']}), flush=True)
if not result['completed']:
    raise RuntimeError('Actual joint path incomplete; last converged accepted pose retained without changing ROM')
