"""Compare exact source-STL vertex connectivity and clustered runtime topology.
Nearest source vertex associations are diagnostics, not a remeshing map."""
import json,os,hashlib,urllib.request
from pathlib import Path
import numpy as np
import trimesh
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
out=Path('test-results/v076-levator-components');out.mkdir(parents=True,exist_ok=True)
data=json.loads(Path('test-results/v076-levator-surface-candidate/levator-binding.json').read_text())
evidence=json.loads(Path('docs/validation/v076/levator-attachment-surfaces.json').read_text())
entry=next(x for x in evidence['sources'] if x['file']=='FMA32540.stl')
cache=Path(os.environ.get('ATLAS_CACHE','/tmp/upper-limb-atlas'));cache.mkdir(parents=True,exist_ok=True)
p=cache/entry['file']
if not p.exists():
    with urllib.request.urlopen('https://raw.githubusercontent.com/Kevin-Mattheus-Moerman/BodyParts3D/2c83c698788cff57f6317bddcf712e7acfbea7d8/assets/BodyParts3D_data/stl/'+entry['file'],timeout=90) as f:p.write_bytes(f.read())
assert hashlib.sha256(p.read_bytes()).hexdigest()==entry['sha256']
raw=trimesh.load_mesh(p,process=False)
v,idx=np.unique(np.asarray(raw.vertices),axis=0,return_inverse=True)
f=idx[np.asarray(raw.faces)]
w=np.asarray(data['rest']).reshape(-1,3);s=data['scale'];o=data['sourceOrigin']
rv=np.column_stack((-w[:,0]/s+o[0],(w[:,2]-.08)/s+o[1],(w[:,1]-1.70)/s+o[2]))
rf=np.asarray(data['faces']).reshape(-1,3)
def audit(v,f):
    edges=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]])
    graph=coo_matrix((np.ones(len(edges)),(edges[:,0],edges[:,1])),shape=(len(v),len(v)))
    n,l=connected_components(graph,directed=False)
    unique,counts=np.unique(np.sort(edges,axis=1),axis=0,return_counts=True)
    rows=[]
    for i in range(n):
        mask=l==i;em=l[unique[:,0]]==i;fm=l[f[:,0]]==i
        assert np.all(l[f[fm]]==i)
        rows.append({'id':i,'vertices':int(mask.sum()),'faces':int(fm.sum()),'boundaryEdges':int(np.sum(em&(counts==1))),'nonmanifoldEdges':int(np.sum(em&(counts>2))),'boundsMm':[v[mask].min(axis=0).tolist(),v[mask].max(axis=0).tolist()]})
    return l,{'vertices':len(v),'faces':len(f),'components':sorted(rows,key=lambda x:-x['vertices'])}
sl,source=audit(v,f);rl,runtime=audit(rv,rf)
dist,nearest=cKDTree(v).query(rv)
for row in runtime['components']:
    mask=rl==row['id'];ids,c=np.unique(sl[nearest[mask]],return_counts=True)
    row['nearestSourceComponents']=[{'id':int(i),'vertices':int(k)} for i,k in zip(ids,c)]
    row['nearestSourceDistanceMaxMm']=float(dist[mask].max())
report={'sourceSha256':entry['sha256'],'sourceExactCoordinateWeld':source,'runtimeClustered':runtime,'productAccepted':False,'limitations':['Vertex-connected components, not measured muscle fascicles','Nearest source vertex correspondence is approximate; not a provenance proof','No deletion, repair, self-intersection check or solid generation']}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print('LEVATOR_COMPONENTS '+json.dumps(report),flush=True)
