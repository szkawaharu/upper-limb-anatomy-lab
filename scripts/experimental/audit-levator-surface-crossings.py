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
def screen(points):
    tri=points[faces];lo=tri.min(1);hi=tri.max(1)
    from rtree.index import Property
    prop=Property();prop.dimension=3
    tree=Index(((i,tuple(np.r_[lo[i],hi[i]]),None) for i in range(len(tri))),properties=prop)
    pairs=[];count=0;parallel_pairs=0;crossings=[]
    def hit(A,B):
        # Moller-Trumbore using full finite edge directions (t in [0,1]).
        e1=B[:,1]-B[:,0];e2=B[:,2]-B[:,0];out=np.zeros(len(A),bool)
        for k in range(3):
            start=A[:,k];d=A[:,(k+1)%3]-start
            h=np.cross(d,e2);det=np.einsum('ij,ij->i',e1,h)
            scale=np.linalg.norm(d,axis=1)*np.linalg.norm(e1,axis=1)*np.linalg.norm(e2,axis=1)
            ok=np.abs(det)>1e-12*np.maximum(scale,1e-30)
            inv=np.divide(1.,det,out=np.zeros_like(det),where=ok)
            s=start-B[:,0];u=np.einsum('ij,ij->i',s,h)*inv
            q=np.cross(s,e1);v=np.einsum('ij,ij->i',d,q)*inv
            t=np.einsum('ij,ij->i',e2,q)*inv
            out|=ok&(u>=-1e-10)&(v>=-1e-10)&(u+v<=1+1e-10)&(t>=-1e-10)&(t<=1+1e-10)
        return out
    def consume(rows):
        nonlocal count,parallel_pairs
        if not rows:return
        ij=np.asarray(rows);A=tri[ij[:,0]];B=tri[ij[:,1]]
        shared=np.any(faces[ij[:,0],:,None]==faces[ij[:,1],None,:],axis=(1,2))
        ij=ij[~shared];A=A[~shared];B=B[~shared];count+=len(ij)
        na=np.cross(A[:,1]-A[:,0],A[:,2]-A[:,0]);nb=np.cross(B[:,1]-B[:,0],B[:,2]-B[:,0])
        parallel_pairs+=int((np.linalg.norm(np.cross(na,nb),axis=1)<=1e-12*np.linalg.norm(na,axis=1)*np.linalg.norm(nb,axis=1)).sum())
        good=hit(A,B)|hit(B,A)
        crossings.extend(map(tuple,ij[good].tolist()))
    for i in range(len(tri)):
        pairs.extend((i,j) for j in tree.intersection(tuple(np.r_[lo[i],hi[i]])) if j>i)
        if len(pairs)>=8192:consume(pairs);pairs=[]
    consume(pairs);tree.close()
    return set(crossings),{'nonadjacentAabbPairs':count,'parallelNormalPairsNotFullyScreened':parallel_pairs}
baseline,br=screen(rest);posed,pr=screen(now)
def shape(p):
    t=p[faces];area=np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)
    edges=np.unique(np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0)
    return area,np.linalg.norm(p[edges[:,1]]-p[edges[:,0]],axis=1)
area0,edge0=shape(rest);area1,edge1=shape(now)
assert area0.min()>0 and edge0.min()>0
report={'index':int(a['index']),'baselineCrossingPairs':len(baseline),'posedCrossingPairs':len(posed),'newCrossingPairs':len(posed-baseline),'newPairs':sorted(posed-baseline),'baselinePairs':sorted(baseline),'baselineScreen':br,'posedScreen':pr,'areaRatioRange':[float((area1/area0).min()),float((area1/area0).max())],'edgeRatioRange':[float((edge1/edge0).min()),float((edge1/edge0).max())],'degenerateFaces':int((area1<=0).sum()),'sharedVertexPairsExcluded':True,'coplanarOverlapCertified':False,'continuousSelfContactCertified':False,'productAccepted':False}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print('LEVATOR_SURFACE_CROSSINGS '+json.dumps(report),flush=True)
