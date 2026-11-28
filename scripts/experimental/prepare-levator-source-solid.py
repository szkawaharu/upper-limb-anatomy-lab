"""Prepare an immutable exact-source levator candidate for offline solid tests.
Only exact index degeneracy and coincident opposite faces may be removed.
All remaining components, including small ones, are recorded and preserved.
The largest closed oriented shell is a separate reference candidate, not a
replacement of the source or a claim of anatomical completeness.
"""
import os,json,hashlib,urllib.request
from pathlib import Path
import numpy as np
import trimesh
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
out=Path('test-results/v076-levator-solid');out.mkdir(parents=True,exist_ok=True)
evidence=json.loads(Path('docs/validation/v076/levator-attachment-surfaces.json').read_text())
entry=next(x for x in evidence['sources'] if x['file']=='FMA32540.stl')
cache=Path(os.environ.get('ATLAS_CACHE','/tmp/upper-limb-atlas'));cache.mkdir(parents=True,exist_ok=True)
path=cache/entry['file']
if not path.exists():
    with urllib.request.urlopen('https://raw.githubusercontent.com/Kevin-Mattheus-Moerman/BodyParts3D/2c83c698788cff57f6317bddcf712e7acfbea7d8/assets/BodyParts3D_data/stl/'+entry['file'],timeout=90) as f:path.write_bytes(f.read())
assert hashlib.sha256(path.read_bytes()).hexdigest()==entry['sha256']
raw=trimesh.load_mesh(path,process=False)
v,idx=np.unique(np.asarray(raw.vertices),axis=0,return_inverse=True)
f=idx[np.asarray(raw.faces)]
np.savez_compressed(out/'original.npz',points=v,faces=f)
groups={};removed=set();pairs=0;degenerate=0
for i,t in enumerate(f):
    if len(set(t))!=3:removed.add(i);degenerate+=1;continue
    groups.setdefault(tuple(sorted(t)),[]).append(i)
def parity(t):return sum(t[i]>t[j] for i in range(3) for j in range(i+1,3))%2
for ids in groups.values():
    plus=[i for i in ids if parity(f[i])==0];minus=[i for i in ids if parity(f[i])==1]
    n=min(len(plus),len(minus));removed.update(plus[:n]);removed.update(minus[:n]);pairs+=n
keep=np.array([i for i in range(len(f)) if i not in removed]);g=f[keep]
used=np.unique(g);remap=np.full(len(v),-1);remap[used]=np.arange(len(used));q=v[used];g=remap[g]
edges=np.concatenate([g[:,[0,1]],g[:,[1,2]],g[:,[2,0]]])
n,labels=connected_components(coo_matrix((np.ones(len(edges)),(edges[:,0],edges[:,1])),shape=(len(q),len(q))),directed=False)
rows=[];shells=[]
for i in range(n):
    vs=np.flatnonzero(labels==i);fs=np.flatnonzero(labels[g[:,0]]==i)
    r=np.full(len(q),-1);r[vs]=np.arange(len(vs));p=q[vs];t=r[g[fs]]
    m=trimesh.Trimesh(p,t,process=False)
    info={'id':i,'vertices':len(p),'faces':len(t),'closed':bool(m.is_watertight),'oriented':bool(m.is_winding_consistent),'signedVolumeMm3':float(m.volume),'minAreaMm2':float(m.area_faces.min()),'euler':int(m.euler_number)}
    rows.append(info);shells.append((p,t,used[vs],keep[fs]))
np.savez_compressed(out/'cleaned-all.npz',points=q,faces=g,originalVertexIds=used,originalFaceIds=keep,removedFaceIds=np.array(sorted(removed)),componentLabels=labels)
main=max(rows,key=lambda x:x['faces']);p,t,vid,fid=shells[main['id']]
report={'sourceSha256':entry['sha256'],'sourceVertices':len(v),'sourceFaces':len(f),'removedOppositePairs':pairs,'removedIndexDegenerateFaces':degenerate,'retainedVertices':len(q),'retainedFaces':len(g),'components':rows,'mainShell':main,'otherShellsPreserved':True,'mainIsExperimentalSubset':True,'sourcePointMovementMm':0,'productAccepted':False}
(out/'preparation.json').write_text(json.dumps(report,indent=2)+'\n')
print('LEVATOR_SOLID_PREP '+json.dumps(report),flush=True)
assert main['closed'] and main['oriented'] and main['signedVolumeMm3']>0 and main['minAreaMm2']>0,'Main shell fails basic geometry gate'
np.savez_compressed(out/'main-shell.npz',points=p,faces=t,originalVertexIds=vid,originalFaceIds=fid)
# Adapt to existing no-boundary-change TetGen driver. Coordinates in atlas mm
# (scale=1 for this independent geometry check); no physical attachment claims.
B={'spec':{'id':'v58Levator'},'rest':p.ravel().tolist(),'faces':t.ravel().tolist(),'weights':[0.5]*len(p),'sourceVertexIds':vid.tolist(),'definition':{'source':{'atlas':'BodyParts3D','worldUnitsPerMm':1,'coordinateSystem':'original atlas mm','file':entry['file']}},'meshPreparation':{'closed':True,'oriented':True,'sourceSha256':entry['sha256'],'mainShellSubset':True,'attachmentWeightsArePlaceholderNotForSolving':True}}
(out/'bindings.json').write_text(json.dumps([B])+'\n')
