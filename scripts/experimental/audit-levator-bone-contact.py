"""Independent surface sample contact audit of saved equilibrium checkpoints.
Compare pointwise penetration depth to the immutable neutral atlas. Vertices,
unique edge midpoints and triangle centroids; not continuous collision proof.
All coordinates in original atlas millimetres, actual browser bone frames."""
import json,os,hashlib,urllib.request
from pathlib import Path
import numpy as np
import trimesh
from scipy.spatial.transform import Rotation
from actual_bone_contact import mesh_distances

src=Path(os.environ['QA_INPUT']);out=Path(os.environ['QA_OUTPUT']);out.mkdir(parents=True,exist_ok=True)
files=sorted(src.glob('accepted-*-index-*.npz'))
assert files
a=np.load(files[0]);faces=a['faces'];n=int(faces.max()+1);rest=a['rest'][:n]+a['center'];classes=a['classes'][:n]
edges=np.unique(np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0)
def samples(v):
    return np.concatenate([v,v[edges].mean(axis=1),v[faces].mean(axis=1)])
fixed=np.concatenate([classes>=0,np.all(classes[edges]>=0,axis=1),np.all(classes[faces]>=0,axis=1)])
p=samples(rest)
frames=json.loads((src/'levator-solid-frames.json').read_text())['frames']
evidence=json.loads(Path('docs/validation/v076/levator-attachment-surfaces.json').read_text())
cache=Path(os.environ.get('ATLAS_CACHE','/tmp/upper-limb-atlas'));cache.mkdir(parents=True,exist_ok=True)
url='https://raw.githubusercontent.com/Kevin-Mattheus-Moerman/BodyParts3D/2c83c698788cff57f6317bddcf712e7acfbea7d8/assets/BodyParts3D_data/stl/'
keys={'C1':'FMA12519.stl','C2':'FMA12520.stl','C3':'FMA12521.stl','C4':'FMA12522.stl','scapula':'FMA13395.stl'}
key=os.environ['QA_BONE'];file=keys[key];path=cache/file
expected=next(x['sha256'] for x in evidence['sources'] if x['file']==file)
if not path.exists():
    with urllib.request.urlopen(url+file,timeout=90) as f:path.write_bytes(f.read())
assert hashlib.sha256(path.read_bytes()).hexdigest()==expected
bone=trimesh.load_mesh(path,process=True)
assert bone.is_watertight,'Bone source not closed'
# Embree accelerates occupancy; nearest-point depth still uses actual triangles.
def query(v):
    ds=[];ins=[]
    for i in range(0,len(v),2048):
        d,_,inside=mesh_distances(bone,v[i:i+2048]);ds.append(d);ins.append(inside)
    return np.concatenate(ds),np.concatenate(ins)
baseline,baseinside=query(p)
Q=np.array([[-1,0,0],[0,0,1],[0,1,0]],float)
rows=[]
for f in files:
    d=np.load(f);index=int(d['index']);assert np.array_equal(d['rest'],a['rest']) and np.array_equal(d['faces'],faces)
    frame=frames[index];s=frame['scale'];o=np.asarray(frame['sourceOrigin']);shift=np.array([0,1.70,.08])
    world=(samples(d['nodes'][:n]+d['center'])-o)@Q.T*s+shift
    b=frame['bones'][key]
    boneworld=Rotation.from_quat(b['q']).inv().apply(world-b['p'])+b['rest']
    current=(boneworld-shift)/s@Q+o
    depth,inside=query(current);increase=depth-baseline;i=int(increase.argmax())
    row={'checkpoint':f.name,'elev':frame['elev'],'neck':frame['neck'],'bone':key,'samples':len(p),'baselineInside':int(baseinside.sum()),'inside':int(inside.sum()),'newInside':int((inside&~baseinside).sum()),'baselineMaxDepthMm':float(baseline.max()),'maxDepthMm':float(depth.max()),'maxIncreaseMm':float(increase[i]),'worstSample':i,'worstKind':'vertex' if i<n else 'edge-midpoint' if i<n+len(edges) else 'triangle-centroid','fixedSampleMaxIncreaseMm':float(increase[fixed].max()),'samplesIncreaseOverPointOneMm':int((increase>.1).sum()),'fixedSamplesIncreaseOverPointOneMm':int((fixed&(increase>.1)).sum())}
    rows.append(row);print('LEVATOR_CONTACT_CASE '+json.dumps(row),flush=True)
report={'bone':key,'boneSHA256':expected,'rows':rows,'sampleCounts':{'vertices':n,'edgeMidpoints':len(edges),'centroids':len(faces)},'engineeringDepthIncreaseGuardMm':.1,'allSamplesWithinGuard':all(r['maxIncreaseMm']<=.1 for r in rows),'productAccepted':False,'continuousContactCertified':False,'selfIntersectionEvaluated':False}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print('LEVATOR_CONTACT_REPORT '+json.dumps(report),flush=True)
