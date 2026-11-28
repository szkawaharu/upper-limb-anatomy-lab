"""Read-only nonadjacent triangle crossing screen on identical source topology.
R-tree AABB candidates, symmetric finite segment/triangle tests. Shared-vertex
pairs and coplanar overlap are excluded; this is NOT a full self-contact proof."""
import os,json
from pathlib import Path
import numpy as np
import trimesh
from rtree.index import Index
src=Path(os.environ['QA_INPUT']);out=Path(os.environ['QA_OUTPUT']);out.mkdir(parents=True,exist_ok=True)
a=np.load(src/os.environ.get('QA_CROSSING_CHECKPOINT','volume-band-candidate.npz'));rest=a['rest'];now=a['nodes'];faces=a['faces']
from surface_crossing_screen import screen
baseline,br=screen(rest,faces);posed,pr=screen(now,faces)
def shape(p):
    t=p[faces];area=np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)
    edges=np.unique(np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0)
    return area,np.linalg.norm(p[edges[:,1]]-p[edges[:,0]],axis=1)
area0,edge0=shape(rest);area1,edge1=shape(now)
assert area0.min()>0 and edge0.min()>0
report={'index':int(a['index']),'baselineCrossingPairs':len(baseline),'posedCrossingPairs':len(posed),'newCrossingPairs':len(posed-baseline),'newPairs':sorted(posed-baseline),'baselinePairs':sorted(baseline),'baselineScreen':br,'posedScreen':pr,'areaRatioRange':[float((area1/area0).min()),float((area1/area0).max())],'edgeRatioRange':[float((edge1/edge0).min()),float((edge1/edge0).max())],'degenerateFaces':int((area1<=0).sum()),'sharedVertexPairsExcluded':True,'coplanarOverlapCertified':False,'continuousSelfContactCertified':False,'productAccepted':False}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print('LEVATOR_SURFACE_CROSSINGS '+json.dumps(report),flush=True)
