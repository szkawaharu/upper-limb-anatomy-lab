"""Experimental attachment candidate on exact exported browser mesh.
Original bones are SHA256-checked. Near-surface masks are not measured footprints.
"""
import os,json,hashlib,urllib.request
from pathlib import Path
import numpy as np
import trimesh
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
out=Path('test-results/v076-levator-surface-candidate');out.mkdir(parents=True,exist_ok=True)
data=json.loads((out/'levator-binding.json').read_text())
evidence=json.loads(Path('docs/validation/v076/levator-attachment-surfaces.json').read_text())
world=np.asarray(data['rest']).reshape(-1,3);s=data['scale'];o=data['sourceOrigin']
p=np.column_stack((-world[:,0]/s+o[0],(world[:,2]-.08)/s+o[1],(world[:,1]-1.70)/s+o[2]))
root=Path(os.environ.get('ATLAS_CACHE','/tmp/upper-limb-atlas'));root.mkdir(parents=True,exist_ok=True)
url='https://raw.githubusercontent.com/Kevin-Mattheus-Moerman/BodyParts3D/2c83c698788cff57f6317bddcf712e7acfbea7d8/assets/BodyParts3D_data/stl/'
def load(file):
    src=next(x for x in evidence['sources'] if x['file']==file)
    path=root/file
    if not path.exists():
        with urllib.request.urlopen(url+file,timeout=90) as f:path.write_bytes(f.read())
    assert hashlib.sha256(path.read_bytes()).hexdigest()==src['sha256'],file
    return trimesh.load_mesh(path,process=True)
files=['FMA12519.stl','FMA12520.stl','FMA12521.stl','FMA12522.stl','FMA13395.stl']
dist=np.stack([np.concatenate([trimesh.proximity.closest_point(m,p[i:i+128])[1] for i in range(0,len(p),128)]) for m in map(load,files)],axis=1)
a=dist[:,:4].min(axis=1);b=dist[:,4];c=dist[:,:4].argmin(axis=1)
f=np.asarray(data['faces']).reshape(-1,3);edges=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]])
_,labels=connected_components(coo_matrix((np.ones(len(edges)),(edges[:,0],edges[:,1])),shape=(len(p),len(p))),directed=False)
summary=[]
# Fixed 1 mm candidate from previous source sensitivity, not a tuned result.
# Compare other radii as coverage diagnostics only.
for radius in [.5,1,2,4]:
    origin=np.flatnonzero((a<=radius)&(a<b));insertion=np.flatnonzero((b<=radius)&(b<a))
    row={'radiusMm':radius,'originCount':len(origin),'insertionCount':len(insertion),'originPerNearestC':np.bincount(c[origin],minlength=4).tolist(),'components':int(labels.max()+1),'componentsWithBoth':len(set(labels[origin])&set(labels[insertion])),'componentsWithNeither':len(set(labels)-set(labels[origin])-set(labels[insertion]))}
    summary.append(row)
    if radius==1:
        result={**row,'rest':data['rest'],'faces':data['faces'],'origin':origin.tolist(),'insertion':insertion.tolist(),'originDistances':a.tolist(),'insertionDistances':b.tolist(),'nearestC':c.tolist(),'sourceSha256':{x['file']:x['sha256'] for x in evidence['sources']},'oldPinCounts':list(map(len,data['pins'])),'productAccepted':False,'anatomicalAcceptance':False}
        assert len(origin) and len(insertion) and all(row['originPerNearestC'])
        (out/'masks.json').write_text(json.dumps(result)+'\n')
report={'runtimeVertices':len(p),'runtimeTriangles':len(f),'oldPinCounts':list(map(len,data['pins'])),'sensitivity':summary,'productAccepted':False,'anatomicalAcceptance':False}
(out/'coverage.json').write_text(json.dumps(report,indent=2)+'\n')
print('LEVATOR_SURFACE_MASKS '+json.dumps(report),flush=True)
