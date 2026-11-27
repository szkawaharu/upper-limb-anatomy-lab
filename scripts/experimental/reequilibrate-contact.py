"""Re-equilibrate an existing pose after audited integration-data changes.
The exact source geometry, fixed masks and final bones must be identical.
No muscle-specific corrective displacement or material coefficient is added.
"""
import sys,json,time,hashlib,shutil
from pathlib import Path
import numpy as np
from actual_bone_contact import ActualBoneContact,atlas_bones
from material_pose import boundary_targets,material_to_world
from bone_sets import model_id,bone_set
from volumetric_reference import Solid
from cage_solver import minimize

root,initial_path=map(Path,sys.argv[1:3]);a=np.load(root/'cage.npz');s=np.load(root/'conforming.npz');old=np.load(initial_path);before=json.loads(initial_path.with_suffix('.json').read_text());clock=time.time()
assert before['progress']==1 and before['allStagesConverged']and before['contactReached']
assert hashlib.sha256(initial_path.read_bytes()).hexdigest()==before['snapshotSha256']
assert np.array_equal(a['nodes'],old['referenceNodes'])and np.array_equal(a['tets'],old['cageTets'])
assert np.array_equal(s['faces'],old['sourceFaces'])and np.array_equal(s['weights'],old['sourceWeights'])
assert np.array_equal(s['embeddingVertices'],old['embeddingVertices'])and np.array_equal(s['embeddingBarycentric'],old['embeddingBarycentric'])
scale=float(a['worldUnitsPerMm']);center=a['centerWorld'];assert np.max(np.abs(s['points']*scale+center-old['referenceWorld']))<1e-12
fixed=np.flatnonzero(s['classes']>=0);solid=Solid(a['nodes'],a['tets'],fixed,element_weights=a['weights']);target=before['target'];initial=old['nodes'];expected=boundary_targets(a['nodes'],fixed,s['classes'],target,center,scale);assert np.max(np.abs(initial[fixed]-expected))<1e-10
identity=model_id(initial_path,before);coverage=bone_set(identity);cal_path=Path('docs/validation/v076/source-bones/browser-numeric/rig-calibration.json');cal=json.loads(cal_path.read_text());contact=ActualBoneContact(s,atlas_bones(coverage,cal,center,scale),center,scale)
code_files=['reequilibrate-contact.py','material_pose.py','actual_bone_contact.py','contact_potential.py','bone_sets.py','contact-bone-sets.json','cage_solver.py','projected_solid_metric.py','constitutive_metric.py','volumetric_reference.py'];code_hashes={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()for name in code_files};(root/'reequilibrium-code-snapshot').mkdir(exist_ok=True)
for name in code_files:shutil.copy2(Path(__file__).with_name(name),root/'reequilibrium-code-snapshot'/name)
def solver(solid,x,**kwargs):
 def observe(nodes,row):
  if row['iteration']%10==0:print(json.dumps({'seconds':time.time()-clock,**row}),flush=True)
 return minimize(solid,x,observer=observe,**kwargs)
x,report=contact.equilibrate(solid,initial,target,solver,bulk=before['bulk'])
points=contact.surface(x);faces=s['faces'];source=s['points'];V=lambda v:float(np.einsum('ti,ti->t',v[faces[:,0]],np.cross(v[faces[:,1]],v[faces[:,2]])).sum()/6)
F=(solid.edges(x)@solid.inverse)[s['elements']];n0=np.cross(source[faces[:,1]]-source[faces[:,0]],source[faces[:,2]]-source[faces[:,0]]);n=np.cross(points[faces[:,1]]-points[faces[:,0]],points[faces[:,2]]-points[faces[:,0]]);expected_normal=np.linalg.det(F)[:,None]*np.einsum('tij,tj->ti',np.linalg.inv(F).transpose(0,2,1),n0)
result={**before,'experimental':True,'accepted':False,'modelId':identity,'contactModelId':identity,'allStagesConverged':report['converged'],'progress':1 if report['converged']else 0,'contactReached':report['converged'],'history':report['history'],'reequilibrium':report,'sourceOpposedFaces':int((np.einsum('ti,ti->t',n,expected_normal)<0).sum()),'surfaceVolumeRatio':V(points)/V(source),'attachmentMaxMm':float(np.linalg.norm(x[fixed]-expected,axis=1).max(initial=0)),'actualDistanceResidualMm':report.get('actualDistanceResidualMm'),'reequilibriumCodeSha256':code_hashes,'inputSha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest()for name in ['cage.npz','conforming.npz','targets.json','preparation.json']},'initialization':{'path':str(initial_path),'snapshotSha256':hashlib.sha256(initial_path.read_bytes()).hexdigest(),'changedData':('unchanged integration control'if hashlib.sha256((root/'cage.npz').read_bytes()).hexdigest()==before['inputSha256']['cage.npz']else 'source/cage material integration weights'),'preserved':['geometry','fixed surface weights','final bones','bulk','stopping tolerance']},'seconds':time.time()-clock,'contactBoneCoverage':coverage,'boneContactEvaluated':True,'selfCollisionEvaluated':False,'activation':None,'force':None}
output=root/initial_path.stem;temporary=output.with_suffix('.tmp.npz');np.savez_compressed(temporary,**{**{k:old[k]for k in old.files},'nodes':x,'sourceWorld':material_to_world(points,target,center,scale),'progress':result['progress']});result['snapshotSha256']=hashlib.sha256(temporary.read_bytes()).hexdigest();temporary.replace(output.with_suffix('.npz'));temporary=output.with_suffix('.tmp.json');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(output.with_suffix('.json'))
print('FINISHED',json.dumps({'converged':report['converged'],'reason':report['reason'],'surfaceVolumeRatio':result['surfaceVolumeRatio'],'sourceOpposedFaces':result['sourceOpposedFaces'],'seconds':time.time()-clock}),flush=True)
if not report['converged']:raise RuntimeError('Re-equilibrium failed; original accepted input retained')
