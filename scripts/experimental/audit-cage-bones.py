"""Actual source bone occupancy/depth, distinct from a contact envelope.
All saved material surface points are tested. Rest overlap is reported separately.
"""
import sys,json,hashlib
from pathlib import Path
import numpy as np,trimesh
from scipy.spatial.transform import Rotation
file=Path(sys.argv[1]);data=np.load(file);report=json.loads(file.with_suffix('.json').read_text())
if report['progress']<1:raise ValueError('Only a completed target may be audited')
cal=json.loads(Path('docs/validation/v076/source-bones/browser-numeric/rig-calibration.json').read_text());scale=cal['scale'];P=np.array([[-1,0,0],[0,0,1],[0,1,0]],float);atlasSC=np.array([-14.91024,-129.66467,1341.37453]);offset=np.array(cal['landmarks']['SC'])/scale-P@atlasSC
rest=data['referenceWorld']/scale;now=data['sourceWorld']/scale;results=[]
for key,name in [('a','FMA13395.stl'),('b','FMA23130.stl')]:
 f=Path('/tmp/upper-limb-atlas')/name;m=trimesh.load_mesh(f,process=True);m.vertices=np.asarray(m.vertices)@P.T+offset
 if not m.is_watertight:raise ValueError('Bone shell not closed: '+name)
 b=report['target'][key];R=Rotation.from_quat(b['q']).as_matrix();query=(now-np.array(b['p'])/scale)@R+np.array(b['rest'])/scale
 before=m.contains(rest);after=m.contains(query);depth_before=np.zeros(len(rest));depth_after=np.zeros(len(rest))
 for points,inside,depth in [(rest,before,depth_before),(query,after,depth_after)]:
  ids=np.flatnonzero(inside)
  for start in range(0,len(ids),512):
   use=ids[start:start+512];depth[use]=trimesh.proximity.closest_point(m,points[use])[1]
 increase=depth_after-depth_before;worst=int(increase.argmax());r={'bone':name,'sourceSha256':hashlib.sha256(f.read_bytes()).hexdigest(),'closed':m.is_watertight,'surfacePoints':len(rest),'restInsidePoints':int(before.sum()),'posedInsidePoints':int(after.sum()),'newInsidePoints':int((after&~before).sum()),'maxRestDepthMm':float(depth_before.max()),'maxPosedDepthMm':float(depth_after.max()),'maxDepthIncreaseMm':float(increase.max()),'worstSurfaceVertex':worst,'increaseOverOneMm':int((increase>1).sum()),'increaseOverPointOneMm':int((increase>.1).sum()),'classification':'vertex-in-shell and vertex-triangle distance; full surface crossings still separate'};results.append(r);print(json.dumps(r),flush=True)
file.with_name(file.stem+'-bones.json').write_text(json.dumps({'experimental':True,'accepted':False,'pose':report['target'],'input':str(file),'inputSha256':hashlib.sha256(file.read_bytes()).hexdigest(),'bones':results},indent=2)+'\n')
