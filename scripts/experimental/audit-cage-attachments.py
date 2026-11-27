"""Report source bone distances for estimated legacy Dirichlet masks.
Distance alone does not identify an anatomical enthesis or validate a mask.
"""
import sys,json,hashlib
from pathlib import Path
import numpy as np,trimesh
p=Path(sys.argv[1]);a=np.load(p/'deltoid-mid-source.npz');calfile=Path('docs/validation/v076/source-bones/browser-numeric/rig-calibration.json');cal=json.loads(calfile.read_text());scale=cal['scale'];P=np.array([[-1,0,0],[0,0,1],[0,1,0]]);offset=np.array(cal['landmarks']['SC'])/scale-P@np.array([-14.91024,-129.66467,1341.37453]);points=a['p']/scale;faces=a['faces'];edges=np.unique(np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0);muscle_spacing=float(np.median(np.linalg.norm(points[edges[:,1]]-points[edges[:,0]],axis=1)));rows=[]
for end,file in [(0,'FMA13395.stl'),(1,'FMA23130.stl')]:
 path=Path('/tmp/upper-limb-atlas')/file;m=trimesh.load_mesh(path,process=True);m.vertices=np.asarray(m.vertices)@P.T+offset;ids=np.flatnonzero(a['weights']==end);near,d,face=trimesh.proximity.closest_point(m,points[ids]);resolution=muscle_spacing+float(np.median(m.edges_unique_length));order=np.argsort(-d);worst=[{'sourceVertex':int(ids[k]),'gapMm':float(d[k]),'nearestBoneFace':int(face[k])}for k in order[:10]]
 row={'end':end,'bone':file,'boneSha256':hashlib.sha256(path.read_bytes()).hexdigest(),'vertices':len(ids),'gapQuantilesMm':np.quantile(d,[0,.1,.5,.9,1]).tolist(),'overTwoMmInvestigationCount':int((d>2).sum()),'measuredResolutionBandMm':resolution,'nearBoneScreeningCount':int((d<=d.min()+resolution).sum()),'screeningIsNotAttachmentClassification':True,'worst':worst}
 if 12473 in ids:
  k=int(np.flatnonzero(ids==12473)[0]);row['collisionExample']={'sourceVertex':12473,'gapMm':float(d[k]),'restPointAxisMm':points[12473].tolist(),'nearestBonePointAxisMm':near[k].tolist(),'nearestBoneFace':int(face[k])}
 rows.append(row);print(json.dumps(row),flush=True)
output={'experimental':True,'accepted':False,'model':'v60DeltMid','method':__doc__,'inputSha256':hashlib.sha256((p/'deltoid-mid-source.npz').read_bytes()).hexdigest(),'calibrationSha256':hashlib.sha256(calfile.read_bytes()).hexdigest(),'muscleMedianEdgeMm':muscle_spacing,'rows':rows};Path(sys.argv[2]).write_text(json.dumps(output,indent=2)+'\n')
