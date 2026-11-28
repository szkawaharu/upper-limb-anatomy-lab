"""Audit linear vertex interpolation between accepted material checkpoints.
Exact cubic determinant extrema (up to floating-point arithmetic), actual
bone-frame attachment error on the .375-degree grid. No contact certificate.
This is a candidate playback screen, not a production interpolation method."""
import os,json
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from volumetric_reference import Solid

def determinant_bounds(F,G):
    D=G-F
    det=lambda a,b,c:np.einsum('ij,ij->i',a,np.cross(b,c))
    a,b,c=F[:,:,0],F[:,:,1],F[:,:,2]
    x,y,z=D[:,:,0],D[:,:,1],D[:,:,2]
    c0=det(a,b,c)
    c1=det(x,b,c)+det(a,y,c)+det(a,b,z)
    c2=det(x,y,c)+det(x,b,z)+det(a,y,z)
    c3=det(x,y,z)
    end=np.linalg.det(G)
    low=np.minimum(c0,end);high=np.maximum(c0,end)
    # Derivative A*t^2+B*t+C. Stable quadratic roots; also handle exact linear case.
    A=3*c3;B=2*c2;C=c1
    disc=B*B-4*A*C
    real=(A!=0)&(disc>=0)
    q=-.5*(B+np.copysign(np.sqrt(np.maximum(disc,0)),B))
    for root in [np.divide(q,A,out=np.full_like(A,np.nan),where=real),
                 np.divide(C,q,out=np.full_like(A,np.nan),where=real&(q!=0)),
                 np.divide(-C,B,out=np.full_like(A,np.nan),where=(A==0)&(B!=0))]:
        active=np.isfinite(root)&(root>0)&(root<1)
        t=root[active];v=((c3[active]*t+c2[active])*t+c1[active])*t+c0[active]
        low[active]=np.minimum(low[active],v);high[active]=np.maximum(high[active],v)
    assert np.isfinite(low).all() and np.isfinite(high).all()
    return low,high

def synthetic():
    # Endpoint J=1; halfway collapse to J=0. A check of endpoints alone misses it.
    f=np.eye(3)[None,:,:];g=np.diag([-1.,-1.,1.])[None,:,:]
    l,h=determinant_bounds(f,g);assert abs(l[0])<1e-12 and abs(h[0]-1)<1e-12
    rng=np.random.default_rng(76)
    f=np.eye(3)[None,:,:]+rng.normal(size=(128,3,3))*.2
    g=np.eye(3)[None,:,:]+rng.normal(size=(128,3,3))*.2
    l,h=determinant_bounds(f,g)
    for t in np.linspace(0,1,101):
        J=np.linalg.det(f+(g-f)*t)
        assert np.all(J>=l-1e-12) and np.all(J<=h+1e-12)
synthetic()
src=Path(os.environ['QA_INPUT']);out=Path(os.environ['QA_OUTPUT']);out.mkdir(parents=True,exist_ok=True)
files=sorted(src.glob('accepted-*-index-*.npz'));assert files
a=np.load(files[0]);rest=a['rest'];fixed=a['fixed'];classes=a['classes'];center=a['center'];tets=a['tets']
solid=Solid(rest,tets,fixed)
frames=json.loads((src/'levator-solid-frames.json').read_text())['frames']
Q=np.array([[-1,0,0],[0,0,1],[0,1,0]],float);keys=['C1','C2','C3','C4','scapula']
def target(i):
    frame=frames[i];s=frame['scale'];o=np.asarray(frame['sourceOrigin']);shift=np.array([0,1.70,.08])
    world=(rest[fixed]+center-o)@Q.T*s+shift;moved=np.empty_like(world)
    for k,key in enumerate(keys):
        mask=classes[fixed]==k;b=frame['bones'][key]
        moved[mask]=Rotation.from_quat(b['q']).apply(world[mask]-b['rest'])+b['p']
    return ((moved-shift)/s)@Q+o-center
old=rest;previous=0;rows=[]
for file in files:
    d=np.load(file)
    for k,v in [('rest',rest),('tets',tets),('fixed',fixed),('classes',classes),('center',center)]:assert np.array_equal(d[k],v)
    index=int(d['index']);new=d['nodes'];assert index!=previous
    F=solid.edges(old)@solid.inverse;G=solid.edges(new)@solid.inverse
    low,high=determinant_bounds(F,G)
    error=0.
    for i in range(min(index,previous),max(index,previous)+1):
        t=(i-previous)/(index-previous)
        error=max(error,float(np.linalg.norm((1-t)*old[fixed]+t*new[fixed]-target(i),axis=1).max()))
    row={'fromIndex':previous,'index':index,'fromDeg':previous*.375,'toDeg':index*.375,'cubicJRange':[float(low.min()),float(high.max())],'elementsOutsideBand':int(((low<.95)|(high>1.05)).sum()),'nonpositiveElements':int((low<=1e-9).sum()),'maxAttachmentErrorMmOnActualFrameGrid':error,'worstMinTet':int(low.argmin()),'worstMaxTet':int(high.argmax())}
    rows.append(row);print('LEVATOR_INTERPOLATION_CASE '+json.dumps(row),flush=True)
    old=new;previous=index
report={'rows':rows,'allSegmentsWithinVolumeGuard':all(r['elementsOutsideBand']==0 for r in rows),'allSegmentsPositive':all(r['nonpositiveElements']==0 for r in rows),'maxAttachmentErrorMm':max(r['maxAttachmentErrorMmOnActualFrameGrid'] for r in rows),'method':'analytic cubic extrema for linear vertex interpolation, sampled actual bone frame attachment error','contactEvaluated':False,'continuousSelfContactCertified':False,'productAccepted':False}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print('LEVATOR_INTERPOLATION_REPORT '+json.dumps(report),flush=True)
