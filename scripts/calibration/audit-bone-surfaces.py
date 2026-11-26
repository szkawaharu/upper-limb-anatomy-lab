"""Independent exact source-triangle proximity; not a cartilage/force solver.
Inside classification uses closed source rib shells. Surface distances are
vertex-to-triangle, never claimed exact triangle-to-triangle gap/contact.
"""
import json,re,sys,hashlib,os
from pathlib import Path
import numpy as np
import trimesh
from scipy.spatial.transform import Rotation
root=Path(sys.argv[1] if len(sys.argv)>1 else '/tmp/upper-limb-atlas')
out=Path(sys.argv[2] if len(sys.argv)>2 else 'test-results/v076-bone-surfaces')
f=json.loads((out/'poses.json').read_text());scale=f['calibration']['scale']
# Work in world-axis millimetres. Translation cancels in distance/classification.
P=np.array([[-1,0,0],[0,0,1],[0,1,0]],float)
offset=np.array(f['calibration']['landmarks']['SC'])/scale-P@f['atlasSC']
manifest=json.loads(re.search(r'const V57_MANIFEST=Object.freeze\((\[.*?\])\);',Path('index.html').read_text()).group(1))
sources=[]
def mesh(name):
    path=root/name;b=path.read_bytes();sources.append({'file':name,'sha256':hashlib.sha256(b).hexdigest()})
    m=trimesh.load_mesh(path,process=True);m.vertices=np.asarray(m.vertices)@P.T+offset;return m
scap=mesh('FMA13395.stl')
ribs=[mesh(x['file'])for x in manifest if x['group']=='ribs'and x['side']=='right']
assert all(m.is_watertight for m in ribs),'Rib sign cannot be interpreted on an open shell'
thorax=trimesh.util.concatenate(ribs)
results=[]
def crossing(mesh,vertices,edges):
    total=0;epsilon=1e-5
    for start in range(0,len(edges),32768):
        e=edges[start:start+32768];a=vertices[e[:,0]];d=vertices[e[:,1]]-a;length=np.linalg.norm(d,axis=1)
        hit,ray,_=mesh.ray.intersects_location(a,d,multiple_hits=False)
        if len(ray):
            t=np.linalg.norm(hit-a[ray],axis=1);total+=int(((t>epsilon)&(t<length[ray]-epsilon)).sum())
    return total
scapEdges=scap.edges_unique;ribEdges=thorax.edges_unique

for pose in f['poses']:
    a=pose['bones']['scapula'];rot=Rotation.from_quat(a['q']).as_matrix();rest=np.array(a['rest'])/scale;at=np.array(a['p'])/scale
    v=(np.asarray(scap.vertices)-rest)@rot.T+at
    # Classify ALL source scapular vertices against the rib union. No mesh
    # simplification, gap threshold relaxation, or proxy shell is used here.
    inside=np.concatenate([thorax.contains(v[i:i+128])for i in range(0,len(v),128)])
    costal=np.array(f['calibration']['contact']['points'])/scale
    costal=(costal-rest)@rot.T+at
    _,distance,_=trimesh.proximity.closest_point(thorax,costal)
    depths=[]
    if inside.any():
        for chunk in np.array_split(v[inside],max(1,int(np.ceil(inside.sum()/128)))):
            depths.extend(trimesh.proximity.closest_point(thorax,chunk)[1].tolist())
    r={'elev':pose['elev'],'plane':pose['plane'],'scapularVertices':len(v),'insideRibVertices':int(inside.sum()),'maxInsideDepthMm':max(depths,default=0),'costalVertexTriangleGapMm':{'min':float(distance.min()),'max':float(distance.max()),'median':float(np.median(distance))},'triangleTriangleIntersections':None}
    if os.environ.get('QA_BONE_EDGES')=='1':
        ribInScap=(np.asarray(thorax.vertices)-at)@rot+rest
        r['surfaceSegmentCrossings']={'scapularEdgesIntoRibs':crossing(thorax,v,scapEdges),'ribEdgesIntoScapula':crossing(scap,ribInScap,ribEdges),'coplanarOverlapEvaluated':False}
    results.append(r);print(json.dumps(r),flush=True)
    (out/'surfaces.json').write_text(json.dumps({'method':__doc__,'sources':sources,'ribShellsClosed':True,'poses':results},indent=2)+'\n')
