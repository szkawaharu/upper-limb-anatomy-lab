"""Actual closed-bone audit at all surface vertices, edge midpoints and face
centroids. Baseline depth is material-point-specific. These finite samples do
not certify all triangle intersections, coplanar overlaps or self-collision.
"""
import sys,json,hashlib
from pathlib import Path
import numpy as np,trimesh
from scipy.spatial.transform import Rotation
from bone_sets import model_id,bone_set
file=Path(sys.argv[1]);a=np.load(file);r=json.loads(file.with_suffix('.json').read_text());assert r['progress']==1
if 'snapshotSha256'in r:assert hashlib.sha256(file.read_bytes()).hexdigest()==r['snapshotSha256']
cal=json.loads(Path('docs/validation/v076/source-bones/browser-numeric/rig-calibration.json').read_text());scale=cal['scale'];P=np.array([[-1,0,0],[0,0,1],[0,1,0]]);offset=np.array(cal['landmarks']['SC'])/scale-P@np.array([-14.91024,-129.66467,1341.37453]);faces=a['sourceFaces'];edges=np.unique(np.sort(np.r_[faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]],axis=1),axis=0);rows=[]
identity=model_id(file,r);coverage=bone_set(identity)
for entry in coverage['bones']:
 key,name=entry['key'],entry['file']
 bonefile=Path('/tmp/upper-limb-atlas')/name;m=trimesh.load_mesh(bonefile,process=True);m.vertices=np.array(m.vertices)@P.T+offset;assert m.is_watertight;bone=r['target'][key];R=Rotation.from_quat(bone['q']).as_matrix();world=a['sourceWorld']/scale;query=(world-np.array(bone['p'])/scale)@R+np.array(bone['rest'])/scale;reference=a['referenceWorld']/scale
 for kind,indices in [('vertex',None),('edge-midpoint',edges),('face-centroid',faces)]:
  rest=reference if indices is None else reference[indices].mean(1);now=query if indices is None else query[indices].mean(1);depths=[];counts=[]
  for points in [rest,now]:
   inside=m.contains(points);ids=np.flatnonzero(inside);depth=np.zeros(len(points))
   for start in range(0,len(ids),256):
    use=ids[start:start+256];depth[use]=trimesh.proximity.closest_point(m,points[use])[1]
   depths.append(depth);counts.append(int(inside.sum()))
  increase=depths[1]-depths[0];worst=int(increase.argmax());row={'bone':name,'samples':kind,'count':len(rest),'insideRest':counts[0],'insidePosed':counts[1],'maxRestDepthMm':float(depths[0].max()),'maxPosedDepthMm':float(depths[1].max()),'maxDepthIncreaseMm':float(increase.max()),'increaseOverOneMm':int((increase>1).sum()),'increaseOverPointOneMm':int((increase>.1).sum()),'worstSample':worst,'worstSourceVertices':[worst]if indices is None else indices[worst].tolist()};rows.append(row);print(json.dumps(row),flush=True)
out={'experimental':True,'accepted':False,'modelId':identity,'boneCoverage':coverage,'input':str(file),'inputSha256':hashlib.sha256(file.read_bytes()).hexdigest(),'method':__doc__,'pose':r['target'],'rows':rows};file.with_name(file.stem+'-surface-samples.json').write_text(json.dumps(out,indent=2)+'\n')
