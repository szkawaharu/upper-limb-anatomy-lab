"""Full repaired source reference regression, NOT production acceptance.
Runs the actual shared solver on immutable BodyParts3D-derived inputs. The
regression motion guard is tied to imposed bone motion and source resolution.
"""
import sys,json,hashlib,tempfile,subprocess,runpy
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'scripts/experimental'))
from snapshot_names import snapshot_name
from volumetric_reference import Solid
from cage_solver import minimize
local_fixture=runpy.run_path(str(root/'tests/reference-local-volume-fixture.py'))
inspect_pose,compare_pressure=local_fixture['inspect_pose'],local_fixture['compare_pressure']
fixture=root/'tests/fixtures/v076-trap-upper-reference';meta=json.loads((fixture/'metadata.json').read_text());file=fixture/'solid.npz';assert hashlib.sha256(file.read_bytes()).hexdigest()==meta['sha256'];a=np.load(file);scale=meta['worldUnitsPerMm'];n=len(a['sourceWeights']);rows=[];positions=[];prescribed=[]
with tempfile.TemporaryDirectory()as temporary:
 p=Path(temporary);src=p/'source';out=p/'reference';src.mkdir();id=meta['id'];np.savez(src/(id+'-source.npz'),p=a['nodesWorld'][:n],faces=a['sourceFaces'],weights=a['sourceWeights']);np.savez(src/(id+'-tetrahedra.npz'),nodes=a['nodesWorld'],tets=a['tets']);(src/'targets.json').write_text(json.dumps({'id':id,'worldUnitsPerMm':scale,'poses':meta['poses']}));subprocess.run([sys.executable,str(root/'scripts/experimental/prepare-source-cage.py'),str(src),str(out)],check=True,capture_output=True)
 c=np.load(out/'cage.npz');fixed=np.flatnonzero(c['classes']>=0);solid=Solid(c['nodes'],c['tets'],fixed,element_weights=c['weights']);local_audits=[];faces=a['sourceFaces'];volume=lambda y:float(np.einsum('ti,ti->t',y[faces[:,0]],np.cross(y[faces[:,1]],y[faces[:,2]])).sum()/6)
 for pose in meta['poses']:
  result=subprocess.run([sys.executable,str(root/'scripts/experimental/run-conforming-reference.py'),str(out),str(pose['elev']),str(pose['plane']),'30'],check=True,capture_output=True,text=True);name=snapshot_name(pose['elev'],pose['plane'],30);r=json.loads((out/(name+'.json')).read_text());s=np.load(out/(name+'.npz'));assert r['progress']==1 and r['allStagesConverged'];assert r['attachmentMaxMm']<=meta['acceptance']['maxAttachmentErrorMm'];assert all(t['cageNonPositive']==0 for t in r['stages']);assert r['sourceOpposedFaces']==0;assert abs(volume(s['sourceWorld'])/volume(s['referenceWorld'])-1)<=meta['acceptance']['maxGlobalVolumeChange'];assert np.isfinite(s['sourceWorld']).all()
  local_audits.append({'elev':pose['elev'],**inspect_pose(solid,s['nodes'])})
  expected=np.zeros((len(fixed),3))
  for key,end in [('a',0),('b',1)]:
   mask=c['classes'][fixed]==end;bone=pose[key];expected[mask]=Rotation.from_quat(bone['q']).apply(a['nodesWorld'][fixed[mask]]-bone['rest'])+bone['p']
  error=float(np.linalg.norm(s['sourceWorld'][fixed]-expected,axis=1).max()/scale);assert error<meta['acceptance']['maxAttachmentErrorMm'];prescribed.append(expected);positions.append(s['sourceWorld']);rows.append({'elev':pose['elev'],'volumeRatio':r['surfaceVolumeRatio'],'attachmentMaxMm':error,'JRange':[r['stages'][-1]['cageJQuantiles'][0],r['stages'][-1]['cageJQuantiles'][-1]],'seconds':r['seconds']})
 maximum=float(np.linalg.norm(positions[1]-positions[0],axis=1).max()/scale);fixed_maximum=float(np.linalg.norm(prescribed[1]-prescribed[0],axis=1).max()/scale);bound=fixed_maximum+meta['acceptance']['sourceMedianEdgeMm'];assert maximum<=bound,'Source deformation exceeds prescribed bone-motion plus measured mesh-resolution guard'
 # Isolated 6-degree material-only comparison; never replace the baseline.
 pressure_experiment=compare_pressure(solid,s['nodes'])
 refinement_experiment=local_fixture['compare_refinement'](solid,s['nodes'],faces)
 mixed_experiment=local_fixture['compare_mixed_pressure'](solid,s['nodes'])
 # Actual iterative return from the loaded 6-degree equilibrium, not reuse
 # of a cached neutral output. Pins return to the immutable atlas-rest data.
 solid=Solid(c['nodes'],c['tets'],fixed,element_weights=c['weights']);initial=solid.predict_boundary(s['nodes'],c['nodes'][fixed]);returned,reset=minimize(solid,initial,bulk=30,tolerance_mm=1e-6,iterations=30);drift=float(np.linalg.norm(returned[:n]-c['nodes'][:n],axis=1).max());assert reset['converged']and drift<.001,'Reference neutral return exceeds 1-micrometre numerical regression guard';assert np.array_equal(returned[fixed],c['nodes'][fixed])
 report={'anatomicalAcceptance':False,'productionEngineAcceptance':False,'fixtureSha256':meta['sha256'],'referenceSourceVertexMaxStepMm':maximum,'prescribedFixedVertexMaxStepMm':fixed_maximum,'engineeringRegressionBoundMm':bound,'referenceNeutralReturn':{'maxDriftMm':drift,'equilibriumResidualMm':reset['residualMaxMm'],'numericalGuardMm':.001,'notAProductionInteractionTest':True},'rows':rows,'localVolumeAudits':local_audits,'localPressureExperiment':pressure_experiment,'interiorRefinementExperiment':refinement_experiment,'mixedPressureExperiment':mixed_experiment,'limitations':meta['limitations']+['Local pressure comparison excludes bone contact and is not anatomical acceptance']};destination=root/'test-results/v076-real-source-reference.json';destination.parent.mkdir(exist_ok=True);destination.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
print('PASS full-source reference: exact real-bone boundary, positive elements, global volume engineering target, no large chart amplification in frozen 3-to-6-degree fixture; not anatomical acceptance')
