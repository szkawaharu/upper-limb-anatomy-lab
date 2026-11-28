"""Source-surface audit for missing levator attachment mask.
Distances are unsigned point-to-triangle distances in the original atlas mm.
No thresholds are anatomical footprints, and no new pins are installed.
"""
import hashlib,json,re,urllib.request,os
from pathlib import Path
import numpy as np
import trimesh
from scipy.sparse.csgraph import connected_components
from scipy.sparse import coo_matrix

out=Path('test-results/v076-levator-attachments');out.mkdir(parents=True,exist_ok=True)
root=Path(os.environ.get('ATLAS_CACHE','/tmp/upper-limb-atlas'));root.mkdir(parents=True,exist_ok=True)
url='https://raw.githubusercontent.com/Kevin-Mattheus-Moerman/BodyParts3D/2c83c698788cff57f6317bddcf712e7acfbea7d8/assets/BodyParts3D_data/stl/'
html=Path('index.html').read_text()
manifest=json.loads(re.search(r'const V57_MANIFEST=Object.freeze\((\[.*?\])\);',html).group(1))
spine=[x for x in manifest if x['id'] in ['C1','C2','C3','C4']]
scap=next(s for s in json.loads(Path('docs/validation/v076/source-bones/after-surfaces-and-edges.json').read_text())['sources'] if s['file']=='FMA13395.stl')
sources=[]
def load(file,sha=None,size=None,sha256=None):
    p=root/file
    if not p.exists():
        with urllib.request.urlopen(url+file,timeout=90) as resp:p.write_bytes(resp.read())
    raw=p.read_bytes();blob=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest();digest=hashlib.sha256(raw).hexdigest()
    if sha:assert sha==blob,(file,'git blob mismatch')
    if sha256:assert sha256==digest,(file,'SHA256 mismatch')
    if size:assert size==len(raw),(file,'size mismatch')
    m=trimesh.load_mesh(p,process=True)
    sources.append({'file':file,'sha256':digest,'gitBlob':blob,'bytes':len(raw),'vertices':len(m.vertices),'faces':len(m.faces)})
    return m
muscle=load('FMA32540.stl',size=1566884)
bones=[load(x['file'],sha=x['sha'],size=x['bytes']) for x in spine]+[load(scap['file'],sha256=scap['sha256'])]
p=np.asarray(muscle.vertices)
def nearest(m):
    return np.concatenate([trimesh.proximity.closest_point(m,p[i:i+128])[1] for i in range(0,len(p),128)])
dist=np.stack([nearest(m) for m in bones],axis=1)
# Exact reconstruction of the two old landmark-segment distances in mm
# on the original unclustered surface. Browser clustering/roundoff is separate.
def segment(a,b):
    a,b=np.asarray(a),np.asarray(b);d=b-a
    t=np.clip((p-a)@d/(d@d),0,1)
    return np.linalg.norm(p-a-t[:,None]*d,axis=1)
old_origin=segment([-33,-66,1480],[-33,-66,1400])
old_insert=segment([-63,-48,1350],[-62,-46,1321])
t=np.clip((old_origin/np.maximum(old_origin+old_insert,1e-9)-.12)/.76,0,1)
w=t*t*(3-2*t)
w[(old_origin<9)&(old_origin<old_insert*.75)]=0
w[(old_insert<8)&(old_insert<old_origin*.75)]=1
def stats(a):
    i=int(np.argmin(a));return {'minMm':float(a[i]),'minVertex':i,'minPointAtlasMm':p[i].tolist(),'percentilesMm':np.percentile(a,[1,5,25,50,95]).tolist()}
edges=muscle.edges_unique
def components(mask):
    ids=np.flatnonzero(mask);e=edges[mask[edges[:,0]]&mask[edges[:,1]]]
    if not len(ids):return []
    inv=np.full(len(p),-1,int);inv[ids]=np.arange(len(ids))
    graph=coo_matrix((np.ones(len(e)),(inv[e[:,0]],inv[e[:,1]])),shape=(len(ids),len(ids)))
    _,labels=connected_components(graph,directed=False)
    return sorted(np.bincount(labels).tolist(),reverse=True)
origin=dist[:,:4].min(axis=1);insert=dist[:,4];nearest_c=dist[:,:4].argmin(axis=1)
sensitivity=[]
for radius in [0.5,1,2,4,8]:
    o=(origin<=radius)&(origin<insert);s=(insert<=radius)&(insert<origin)
    sensitivity.append({'radiusMm':radius,'originVertices':int(o.sum()),'insertionVertices':int(s.sum()),'originPerNearestC':[int((o&(nearest_c==i)).sum()) for i in range(4)],'originComponents':components(o),'insertionComponents':components(s),'overlapNearBoth':int(((origin<=radius)&(insert<=radius)).sum()),'originBoundsAtlasMm':[p[o].min(axis=0).tolist(),p[o].max(axis=0).tolist()] if o.any() else None,'insertionBoundsAtlasMm':[p[s].min(axis=0).tolist(),p[s].max(axis=0).tolist()] if s.any() else None})
report={'sources':sources,'sourceMuscleBoundsAtlasMm':muscle.bounds.tolist(),'legacyLandmarks':{'origin':stats(old_origin),'insertion':stats(old_insert),'fixedOrigin':int((w==0).sum()),'fixedInsertion':int((w==1).sum()),'weightRange':[float(w.min()),float(w.max())]},'sourceBoneDistances':{k:stats(dist[:,i]) for i,k in enumerate(['C1','C2','C3','C4','scapula'])},'nearSurfaceSensitivity':sensitivity,'productAccepted':False,'anatomicalAcceptance':False,'limitations':['Unsigned whole-bone distance, not attachment footprint or transverse-process segmentation','Unclustered source surface, not browser binding','No new pins selected']}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print('LEVATOR_ATTACHMENT_SURFACES '+json.dumps(report),flush=True)
