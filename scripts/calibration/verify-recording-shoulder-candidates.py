"""Independent source geometry check of limited SC/AC candidate poses.
No global convex hull, no new physiological gap targets and no product switch.
All original scapular vertices are tested per closed rib shell, with both
edge-to-face directions. Coplanar contact is not classified.
"""
import hashlib,json,os,sys,urllib.request
from pathlib import Path
import numpy as np
import trimesh
from scipy.spatial.transform import Rotation

out=Path('test-results/v076-recording-shoulder')
data=json.loads((out/'candidates.json').read_text())
cal=data['calibration'];scale=cal['scale']
archive=json.loads(Path('docs/validation/v076/source-bones/after-surfaces-and-edges.json').read_text())
sources=archive['sources']
root=Path(os.environ.get('ATLAS_CACHE','/tmp/upper-limb-atlas'));root.mkdir(parents=True,exist_ok=True)
url='https://raw.githubusercontent.com/Kevin-Mattheus-Moerman/BodyParts3D/2c83c698788cff57f6317bddcf712e7acfbea7d8/assets/BodyParts3D_data/stl/'
P=np.array([[-1,0,0],[0,0,1],[0,1,0]],float)
offset=np.array(cal['landmarks']['SC'])/scale-P@data['atlasSC']
def load(source):
    path=root/source['file']
    if not path.exists():
        with urllib.request.urlopen(url+source['file'],timeout=90) as resp:
            path.write_bytes(resp.read())
    if hashlib.sha256(path.read_bytes()).hexdigest()!=source['sha256']:
        raise ValueError('Source hash mismatch: '+source['file'])
    m=trimesh.load_mesh(path,process=True)
    m.vertices=np.asarray(m.vertices)@P.T+offset
    return m
scap=load(next(s for s in sources if s['file']=='FMA13395.stl'))
ribs=[load(s) for s in sources if s['file']!='FMA13395.stl']
assert len(ribs)==12 and all(m.is_watertight for m in ribs)
thorax=trimesh.util.concatenate(ribs)
assert scap.is_watertight,'Source scapula sign classification requires a closed shell'
print(json.dumps({'sourceScapulaVertices':len(scap.vertices),'ribs':len(ribs),'acceleratedRay':bool(trimesh.ray.has_embree)}),flush=True)
assert len(scap.vertices)==56630
costal=np.asarray(cal['contact']['points'],float)/scale
assert costal.ndim==2 and costal.shape[1]==3
epsilon=1e-5
def inside(mesh,v):
    lo,hi=mesh.bounds
    ids=np.flatnonzero(np.all((v>=lo)&(v<=hi),axis=1))
    hit=np.zeros(len(v),bool)
    for chunk in np.array_split(ids,max(1,int(np.ceil(len(ids)/512)))):
        if len(chunk):hit[chunk]=mesh.contains(v[chunk])
    return hit
def crossings(mesh,v,edges):
    lo,hi=mesh.bounds
    n=0
    for start in range(0,len(edges),8192):
        e=edges[start:start+8192];a=v[e[:,0]];b=v[e[:,1]]
        mask=np.all((np.minimum(a,b)<=hi)&(np.maximum(a,b)>=lo),axis=1)
        a,b=a[mask],b[mask]
        if not len(a):continue
        d=b-a;l=np.linalg.norm(d,axis=1)
        good=l>2*epsilon;a,d,l=a[good],d[good],l[good]
        if not len(a):continue
        hit,ray,_=mesh.ray.intersects_location(a,d,multiple_hits=True)
        if len(ray):
            t=np.linalg.norm(hit-a[ray],axis=1)
            n+=int(np.unique(ray[(t>epsilon)&(t<l[ray]-epsilon)]).size)
    return n
rest_gap=None;rows=[]
for row in data['rows']:
    a=row['bones']['scapula'];rot=Rotation.from_quat(a['q']).as_matrix()
    rest=np.array(a['rest'])/scale;at=np.array(a['p'])/scale
    v=(np.asarray(scap.vertices)-rest)@rot.T+at
    s=(costal-rest)@rot.T+at
    _,dist,_=trimesh.proximity.closest_point(thorax,s)
    if rest_gap is None:rest_gap=dist.copy()
    # Complete source vertex classification, not just the 987 representatives.
    hits=np.zeros(len(v),bool)
    rib_in_scap=0
    for rib in ribs:
        hits|=inside(rib,v)
        inverse=(np.asarray(rib.vertices)-at)@rot+rest
        rib_in_scap+=int(inside(scap,inverse).sum())
    # Both directions on the original triangle edges. If vertex penetration
    # already rejects this candidate we still record it and do not claim
    # the unevaluated edge count as zero.
    if hits.any() or rib_in_scap:
        cross=None
    else:
        cross={'scapularEdgesIntoRibs':0,'ribEdgesIntoScapula':0}
        for rib in ribs:
            inverse=(np.asarray(rib.vertices)-at)@rot+rest
            cross['scapularEdgesIntoRibs']+=crossings(rib,v,scap.edges_unique)
            cross['ribEdgesIntoScapula']+=crossings(scap,inverse,rib.edges_unique)
    r={**{k:row[k] for k in ['elev','plane','alpha','parameters','priorDeviationL2Deg','controlBackend']},
       'sourceScapularVertices':len(v),'insideRibVertices':int(hits.sum()),'insideScapulaRibVertices':rib_in_scap,
       'surfaceSegmentCrossings':cross,'coplanarOverlapEvaluated':False,
       'costalVertexTriangleGapMm':{'min':float(dist.min()),'median':float(np.median(dist)),'max':float(dist.max()),
          'rmsDeltaFromRest':float(np.sqrt(np.mean((dist-rest_gap)**2)))},
       'screenedCollisionFree':bool(not hits.any() and not rib_in_scap and cross is not None and not any(cross.values())),
       'productAccepted':False}
    rows.append(r);print(json.dumps(r),flush=True)
    (out/'report.json').write_text(json.dumps({'sourceSha256':sources,'rows':rows,'anatomicalAcceptance':False,'productAccepted':False},indent=2)+'\n')
# Within each posture, keep only collision-screened candidates with smaller
# rest-relative costal RMS than the unchanged conservative control.
selected=[]
for e,p in sorted({(r['elev'],r['plane']) for r in rows}):
    group=[r for r in rows if r['elev']==e and r['plane']==p];control=group[0]
    feasible=[r for r in group[1:] if r['screenedCollisionFree'] and r['costalVertexTriangleGapMm']['rmsDeltaFromRest']<control['costalVertexTriangleGapMm']['rmsDeltaFromRest'] and r['priorDeviationL2Deg']<=control['priorDeviationL2Deg']]
    selected.append({'elev':e,'plane':p,'control':control,'candidate':min(feasible,key=lambda r:r['costalVertexTriangleGapMm']['rmsDeltaFromRest']) if feasible else None})
report={'sourceSha256':sources,'rows':rows,'selected':selected,'anatomicalAcceptance':False,'productAccepted':False,
        'limitations':['Limited poses and line segment of six linked parameters only','Vertex-triangle distances, not soft tissue measurements','No coplanar classification, cartilage, muscles or continuous motion acceptance']}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'selected':selected,'productAccepted':False}),flush=True)
