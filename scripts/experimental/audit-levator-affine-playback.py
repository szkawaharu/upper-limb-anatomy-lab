"""Offline reduced affine boundary correction for playback.
Five rigid bone transforms each supply 12 affine coefficients. A precomputed
linear coupled FEM extension propagates their residual from interpolated
checkpoint transforms. No nonlinear solve at playback. Sampled screening,
not equilibrium, continuous collision certification or a released renderer.
"""
import os,json,time
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.sparse.linalg import splu
from volumetric_reference import Solid
from volumetric_preconditioner import strain_metric_matrix

src=Path(os.environ['QA_INPUT']);out=Path(os.environ['QA_OUTPUT']);out.mkdir(parents=True,exist_ok=True)
files=sorted(src.glob('accepted-*-index-*.npz'));assert len(files)==20
a=np.load(files[0]);rest=a['rest'];fixed=a['fixed'];classes=a['classes'];center=a['center'];tets=a['tets'];faces=a['faces']
solid=Solid(rest,tets,fixed);assert len(solid.inactive)==0
frames=json.loads((src/'levator-solid-frames.json').read_text())['frames']
assert len(frames)==81 and all(f['index']==i and f['elev']==i*.375 for i,f in enumerate(frames))
Q=np.array([[-1,0,0],[0,0,1],[0,1,0]],float);keys=['C1','C2','C3','C4','scapula']
scale=100.
features=np.c_[np.ones(len(rest)),rest/scale]
probe=np.vstack([np.zeros(3),np.eye(3)*scale])
def transform(points,frame,key):
    s=frame['scale'];o=np.asarray(frame['sourceOrigin']);shift=np.array([0,1.70,.08])
    world=(points+center-o)@Q.T*s+shift;b=frame['bones'][key]
    moved=Rotation.from_quat(b['q']).apply(world-b['rest'])+b['p']
    return ((moved-shift)/s)@Q+o-center
def coefficients(frame):
    result=[]
    for key in keys:
        v=transform(probe,frame,key)-probe
        result.append(np.vstack([v[0],v[1:]-v[0]]))
    return np.asarray(result).ravel()
coeff=np.asarray([coefficients(f) for f in frames])
# 60 full affine modes, including translations, in centered mm coordinates.
P=np.zeros((len(fixed),3,60))
for i in range(5):
    mask=classes[fixed]==i
    for k in range(4):
        for c in range(3):
            P[mask,c,i*12+k*3+c]=features[fixed[mask],k]
def target(i):
    d=rest[fixed]+(P@coeff[i])
    actual=np.empty_like(d)
    for k,key in enumerate(keys):
        mask=classes[fixed]==k
        actual[mask]=transform(rest[fixed[mask]],frames[i],key)
    assert np.max(np.linalg.norm(actual-d,axis=1))<1e-9
    return actual
targets=[target(i) for i in range(len(frames))]
start=time.monotonic()
free=(solid.free[:,None]*3+np.arange(3)).ravel();pins=(fixed[:,None]*3+np.arange(3)).ravel()
H=strain_metric_matrix(solid,rest,1000,full=True)
A=H[free][:,free].tocsc();rhs=-(H[free][:,pins]@P.reshape(-1,60))
factor=splu(A,permc_spec='MMD_AT_PLUS_A',diag_pivot_thresh=0)
X=factor.solve(rhs)
error=float(np.max(np.abs(A@X-rhs))/max(float(np.asarray(abs(A).sum(axis=1)).max())*np.max(np.abs(X))+np.max(np.abs(rhs)),1e-300))
assert np.isfinite(X).all() and error<1e-10
basis=np.zeros((len(rest)*3,60));basis[free]=X;basis[pins]=P.reshape(-1,60)
# Explicit test of common translation: constant translation across all five labels.
for c in range(3):
    v=np.zeros(60);v[c::12]=1.
    expected=np.zeros_like(rest);expected[:,c]=1.
    assert np.max(np.abs((basis@v).reshape(-1,3)-expected))<1e-7
basis32=basis.astype(np.float32)
np.savez_compressed(out/'affine-correction-basis.npz',basis=basis32.reshape(len(rest),3,60),featuresScaleMm=scale,keys=keys,fixed=fixed,classes=classes,center=center)
print('LEVATOR_AFFINE_SETUP '+json.dumps({'modes':60,'bytesFloat32':basis32.nbytes,'linearBackwardError':error,'seconds':time.monotonic()-start}),flush=True)
def evaluate(z,i):
    J=np.linalg.det(solid.edges(z)@solid.inverse)
    return {'minJ':float(J.min()),'maxJ':float(J.max()),'outside':int(((J<.95)|(J>1.05)).sum()),'nonPositive':int((J<=1e-9).sum()),'pinErrorMm':float(np.linalg.norm(z[fixed]-targets[i],axis=1).max())}
old=rest;previous=0;rows=[]
for file in files:
    d=np.load(file)
    for k,v in [('rest',rest),('tets',tets),('fixed',fixed),('classes',classes),('center',center),('faces',faces)]:assert np.array_equal(d[k],v)
    index=int(d['index']);new=d['nodes'];assert index!=previous
    for i in range(min(index,previous),max(index,previous)+1):
        t=(i-previous)/(index-previous);linear=(1-t)*old+t*new
        delta=coeff[i]-((1-t)*coeff[previous]+t*coeff[index])
        corrected=linear+(basis@delta).reshape(-1,3)
        # Simulate float32 basis and coefficient storage, then exact pin overwrite.
        z32=linear+(basis32@delta.astype(np.float32)).reshape(-1,3)
        z32[fixed]=targets[i]
        row={'fromIndex':previous,'toIndex':index,'index':i,'linear':evaluate(linear,i),'corrected':evaluate(corrected,i),'float32BasisExactPins':evaluate(z32,i),'maxBasisQuantizationErrorMm':float(np.linalg.norm(z32-corrected,axis=1).max())}
        rows.append(row)
    old=new;previous=index
def summarize(key):
    r=[v[key] for v in rows]
    return {'minJ':min(v['minJ'] for v in r),'maxJ':max(v['maxJ'] for v in r),'samplePosesOutsideBand':sum(v['outside']>0 for v in r),'maxOutsideElements':max(v['outside'] for v in r),'maxNonPositiveElements':max(v['nonPositive'] for v in r),'maxPinErrorMm':max(v['pinErrorMm'] for v in r)}
report={'sampleCount':len(rows),'modes':60,'basisFloat32Bytes':basis32.nbytes,'linearBackwardError':error,'linear':summarize('linear'),'corrected':summarize('corrected'),'float32BasisExactPins':summarize('float32BasisExactPins'),'rows':rows,'continuousJCertified':False,'contactEvaluated':False,'equilibriumEvaluated':False,'productAccepted':False}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print('LEVATOR_AFFINE_REPORT '+json.dumps({k:v for k,v in report.items() if k!='rows'}),flush=True)
