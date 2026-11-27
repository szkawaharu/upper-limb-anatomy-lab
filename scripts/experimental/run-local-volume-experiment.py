"""Opt-in 90-degree local-volume experiment using an existing audited cage.
No writes to the input directory or accepted pose. Failed solve writes a report
only. Even a converged experiment is not anatomical/production acceptance.
"""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from actual_bone_contact import ActualBoneContact, atlas_bones
from material_pose import boundary_targets, material_to_world
from bone_sets import model_id, bone_set
from volumetric_reference import Solid
from local_volume_constraints import equilibrate, physical_cells, audit

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cage_directory', type=Path)
    parser.add_argument('snapshot', type=Path)
    parser.add_argument('output_directory', type=Path)
    parser.add_argument('--max-log-error', required=True, type=float,
                        help='Experimental maximum absolute log(J), not release tolerance')
    parser.add_argument('--rho', type=float, default=30.)
    parser.add_argument('--outer-steps', type=int, default=20)
    parser.add_argument('--audit-only', action='store_true')
    args = parser.parse_args()
    root, initial_path, out = args.cage_directory, args.snapshot, args.output_directory
    if out.resolve() == root.resolve() or out.exists():
        raise ValueError('Use a new output directory; existing results never overwritten')
    if not np.isfinite(args.max_log_error) or args.max_log_error <= 0:
        raise ValueError('Invalid local experimental tolerance')
    a = np.load(root/'cage.npz')
    s = np.load(root/'conforming.npz')
    old = np.load(initial_path)
    report_path = initial_path.with_suffix('.json')
    previous = json.loads(report_path.read_text())
    assert previous['progress'] == 1 and previous['allStagesConverged'] and previous['contactReached']
    assert previous['snapshotSha256'] == sha(initial_path)
    assert np.array_equal(a['nodes'], old['referenceNodes']) and np.array_equal(a['tets'], old['cageTets'])
    for key, other in [('faces','sourceFaces'),('weights','sourceWeights'),
                       ('embeddingVertices','embeddingVertices'),('embeddingBarycentric','embeddingBarycentric')]:
        assert np.array_equal(s[key],old[other]), key
    assert np.array_equal(a['classes'], s['classes'])
    scale, center = float(a['worldUnitsPerMm']), a['centerWorld']
    assert np.max(np.abs(s['points']*scale+center-old['referenceWorld'])) < 1e-12
    fixed = np.flatnonzero(s['classes'] >= 0)
    solid = Solid(a['nodes'], a['tets'], fixed, a['weights'])
    target, initial = previous['target'], old['nodes']
    expected = boundary_targets(a['nodes'], fixed, s['classes'], target, center, scale)
    assert np.max(np.abs(initial[fixed]-expected), initial=0) < 1e-10
    cells,_ = physical_cells(solid,a['sourceMaterialVolumes'],a['backgroundWeights'])
    cal_path = Path('docs/validation/v076/source-bones/browser-numeric/rig-calibration.json')
    inputs = [root/'cage.npz',root/'conforming.npz',root/'preparation.json',initial_path,report_path,cal_path]
    code = ['run-local-volume-experiment.py','local_volume_constraints.py','actual_bone_contact.py',
            'material_pose.py','bone_sets.py','contact-bone-sets.json','volumetric_reference.py',
            'cage_solver.py','projected_solid_metric.py','constitutive_metric.py','contact_potential.py']
    result = {'experimental':True,'accepted':False,'target':target,'initialAudit':audit(solid,initial,cells),
              'inputsSha256':{str(p):sha(p) for p in inputs},
              'codeSha256':{p:sha(Path(__file__).with_name(p)) for p in code},
              'parameters':{'maxAbsLogJ':args.max_log_error,'rho':args.rho,'outerSteps':args.outer_steps,
                            'baseBulk':previous['bulk']},
              'productionEngineAcceptance':False,'anatomicalAcceptance':False,'selfCollisionEvaluated':False}
    out.mkdir(parents=True,exist_ok=False)
    if not args.audit_only:
        coverage = bone_set(model_id(initial_path,previous))
        bones = atlas_bones(coverage,json.loads(cal_path.read_text()),center,scale)
        result['boneInputs'] = [{k:b[k] for k in ('key','file','sha256')} for b in bones]
        contact = ActualBoneContact(s,bones,center,scale)
        x, solve = equilibrate(solid,initial,contact,target,a['sourceMaterialVolumes'],a['backgroundWeights'],
                               max_log_error=args.max_log_error,rho=args.rho,bulk=previous['bulk'],
                               outer_steps=args.outer_steps)
        result['solve'] = solve
        if solve['converged']:
            points=contact.surface(x)
            result['independentBoneScan']=contact.scan(x,target)[0]
            result['attachmentMaxMm']=float(np.linalg.norm(x[fixed]-expected,axis=1).max(initial=0))
            result['finalAudit']=audit(solid,x,cells)
            path=out/'experimental-candidate.npz'
            # Intentionally not an accepted joint-path snapshot schema.
            np.savez_compressed(path,nodes=x,sourceWorld=material_to_world(points,target,center,scale),
                                referenceNodes=a['nodes'],cageTets=a['tets'])
            result['candidateSha256']=sha(path)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__ == '__main__':
    main()
