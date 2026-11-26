"""Locate missing fixed ends and source-surface distances. Minimum-distance
regions are candidates, not histologically measured anatomical footprints.
"""
import json,re
from pathlib import Path
import numpy as np,trimesh
src=Path('test-results/v076-source-clean/bindings.json');root=Path('/tmp/upper-limb-atlas');fixture=json.load(open('docs/validation/v076/contact-rms/rig-calibration.json'));scale=fixture['scale'];SC=np.array([-14.91024,-129.66467,1341.37453]);P=np.array([[-1,0,0],[0,0,1],[0,1,0]],float);offset=np.array(fixture['landmarks']['SC'])/scale-P@SC
bones={k:trimesh.load_mesh(root/f,process=True)for k,f in {'scapula':'FMA13395.stl','clavicle':'FMA13322.stl','humerus':'FMA23130.stl'}.items()}
records=[]
for B in json.loads(src.read_text()):
 p=(np.array(B['rest']).reshape(-1,3)/scale-offset)@P;w=np.array(B['weights']);r={'id':B['spec']['id'],'name':B['spec']['name'],'fixedOrigin':int((w==0).sum()),'fixedInsertion':int((w==1).sum()),'weightRange':[float(w.min()),float(w.max())],'ends':[]}
 for end,bone in enumerate([B['spec']['a'],B['spec']['b']]):
  ids=np.argsort(w if end==0 else-w)[:max(1,int(len(w)*.01))];points=p[ids];entry={'bone':bone,'centroidMm':points.mean(0).tolist(),'candidateBoxMm':[points.min(0).tolist(),points.max(0).tolist()],'regionMethod':'1% extreme weights for audit only; not a new attachment mask'}
  if bone in bones:
   closest,d,_=trimesh.proximity.closest_point(bones[bone],points);entry['boneVertexTriangleDistanceMm']={'min':float(d.min()),'max':float(d.max()),'median':float(np.median(d))};entry['nearestBonePatchMm']=closest.mean(0).tolist()
  r['ends'].append(entry)
 records.append(r);print(json.dumps(r),flush=True)
Path('test-results/v076-source-clean/attachment-audit.json').write_text(json.dumps({'method':__doc__,'muscles':records},indent=2)+'\n')
