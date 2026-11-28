"""Frozen-pose comparison of stretch regularization with the local J band.
Fixed bone positions and reference integration unchanged. No anatomical
acceptance inferred from an experimental shape regularizer."""
import os,json,time
from pathlib import Path
import numpy as np
from volumetric_reference import Solid
from local_volume_band import VolumeBand
from local_volume_constraints import log_volume_jacobian
from local_stretch_band import StretchBand
from contact_potential import AugmentedSolid
from cage_solver import minimize

src=Path(os.environ['QA_INPUT']);out=Path(os.environ['QA_OUTPUT']);out.mkdir(parents=True,exist_ok=True)
index=int(os.environ['QA_INDEX']);files=list(src.glob(f'accepted-*-index-{index:03d}.npz'));assert len(files)==1
a=np.load(files[0]);rest=a['rest'];tets=a['tets'];fixed=a['fixed'];faces=a['faces'];x=a['nodes'].copy()
solid=Solid(rest,tets,fixed);pins=x[fixed].copy();bulk=1000.
initial_multipliers=a['multipliers'].copy();assert initial_multipliers.shape==(2,len(tets))
edges=np.unique(np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0)
def area(v):
    p=v[faces];return np.linalg.norm(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]),axis=1)
def length(v):return np.linalg.norm(v[edges[:,1]]-v[edges[:,0]],axis=1)
area0=area(rest);edge0=length(rest)
def metrics(v):
    E,g,J=solid.evaluate(v,bulk);assert np.isfinite(E) and J.min()>1e-9
    s=np.linalg.svd(solid.edges(v)@solid.inverse,compute_uv=False)
    ar=area(v)/area0;er=length(v)/edge0
    return {'JRange':[float(J.min()),float(J.max())],'outsideVolumeCount':int(((J<.95)|(J>1.05)).sum()),'principalStretchRange':[float(s.min()),float(s.max())],'maxConditionF':float(np.max(s[:,0]/s[:,-1])),'stretchOutsideCount':int(((s<.5)|(s>2)).any(axis=1).sum()),'areaRatioRange':[float(ar.min()),float(ar.max())],'edgeRatioRange':[float(er.min()),float(er.max())],'attachmentErrorMm':float(np.linalg.norm(v[fixed]-pins,axis=1).max()),'baseEnergy':E}
initial=metrics(x)
assert initial['outsideVolumeCount']==0
band=VolumeBand(solid,300)
_,A=log_volume_jacobian(solid,x,band.cells);_,g,_=solid.evaluate(x,bulk)
g+=(A.T@(band.weights*(initial_multipliers[1]-initial_multipliers[0]))).reshape(x.shape)
assert np.linalg.norm(solid.factor.solve(g[solid.free]),axis=1).max()<.001
rows=[]
p=initial_multipliers
for strength in [.1,1.,10.,100.]:
    z=x.copy();multipliers=p.copy();stages=[]
    stretch=StretchBand(solid,strength)
    for rho in [300,1000,3000,10000,10000,10000]:
        start=time.monotonic();band=VolumeBand(solid,rho,multipliers=multipliers)
        augmented=AugmentedSolid(solid,[band,stretch])
        trial,info=minimize(augmented,z,bulk=bulk,iterations=60,tolerance_mm=.001)
        m=metrics(trial);assert np.array_equal(trial[fixed],pins)
        _,g,_=augmented.evaluate(trial,bulk)
        residual=float(np.linalg.norm(solid.factor.solve(g[solid.free]),axis=1).max())
        row={'strength':strength,'rho':rho,'converged':bool(info['converged'] and residual<.001),'residualMm':residual,'iterations':info['iterations'],'seconds':time.monotonic()-start,**m}
        stages.append(row);print('LEVATOR_STRETCH_CASE '+json.dumps(row),flush=True)
        if not row['converged']:break
        z=trial;multipliers=band.update(z)
        if m['outsideVolumeCount']==0:
            x=z;p=multipliers
            folder=out/f'strength-{strength:g}';folder.mkdir(exist_ok=True)
            np.savez_compressed(folder/'volume-band-candidate.npz',nodes=x,rest=rest,tets=tets,faces=faces,fixed=fixed,classes=a['classes'],index=index,center=a['center'],multipliers=p)
            break
    ok=stages[-1]['converged'] and stages[-1]['outsideVolumeCount']==0
    rows.append({'strength':strength,'candidateSaved':ok,'stages':stages})
    if not ok:break
report={'index':index,'initial':initial,'rows':rows,'stretchBandDiagnostic':[.5,2.],'productAccepted':False,'crossingsEvaluated':False,'contactEvaluated':False}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print('LEVATOR_STRETCH_REPORT '+json.dumps(report),flush=True)
if not all(r['candidateSaved'] for r in rows):raise SystemExit(2)
