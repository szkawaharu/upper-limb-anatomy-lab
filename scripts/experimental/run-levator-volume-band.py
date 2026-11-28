"""Frozen 15-degree local volume-band experiment, no production adoption."""
import os,json,time
from pathlib import Path
import numpy as np
from volumetric_reference import Solid
from cage_solver import minimize
from contact_potential import AugmentedSolid
from local_volume_band import VolumeBand
from local_volume_constraints import audit

src=Path(os.environ['QA_INPUT']);out=Path(os.environ['QA_OUTPUT']);out.mkdir(parents=True,exist_ok=True)
a=np.load(src/os.environ.get('QA_BAND_CHECKPOINT','bulk-3000.npz'));index=int(a['index']);assert 0<=index<=40
rest=a['rest'];tets=a['tets'];fixed=a['fixed'];x=a['nodes'].copy()
solid=Solid(rest,tets,fixed);pins=x[fixed].copy()
bulk=float(os.environ.get('QA_BAND_BULK','3000'));faces=a['faces'];assert np.isfinite(bulk) and bulk>0
def metrics(z):
    E,g,J=solid.evaluate(z,bulk)
    assert np.isfinite(E) and J.min()>1e-9
    p=z[faces];q=rest[faces]
    area=lambda v:np.linalg.norm(np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0]),axis=1)
    ratio=area(p)/area(q)
    ids=np.flatnonzero((J<.95)|(J>1.05))
    return {'JRange':[float(J.min()),float(J.max())],'outsideCount':len(ids),'outsideVolumeFraction':float(solid.weights[ids].sum()),'volumeRatio':float(solid.weights@J),'surfaceAreaRatioRange':[float(ratio.min()),float(ratio.max())],'baseEnergy':E,'baseResidualMm':float(np.linalg.norm(solid.factor.solve(g[solid.free]),axis=1).max()),'attachmentErrorMm':float(np.linalg.norm(z[fixed]-pins,axis=1).max())}
initial=metrics(x);assert initial['baseResidualMm']<.001
diagnosis=audit(solid,x,np.arange(len(tets)),max_iterations=512)
diagnosis['immovableCellCount']=len(diagnosis.pop('immovableCellIds'))
print('LEVATOR_BAND_FEASIBILITY '+json.dumps(diagnosis),flush=True)
rows=[];multipliers=None;accepted=initial['outsideCount']==0
final_rho=None
for k,rho in enumerate([30,100,300,1000,3000,10000]):
    if accepted:break
    start=time.monotonic()
    band=VolumeBand(solid,rho,multipliers=multipliers)
    augmented=AugmentedSolid(solid,[band])
    z,info=minimize(augmented,x,bulk=bulk,iterations=40,tolerance_mm=.001)
    row={'rho':rho,'converged':info['converged'],'reason':info['reason'],'iterations':info['iterations'],'augmentedResidualMm':info['residualMaxMm'],**metrics(z),'seconds':time.monotonic()-start}
    assert np.array_equal(z[fixed],pins)
    rows.append(row);print('LEVATOR_BAND_CASE '+json.dumps(row),flush=True)
    if not info['converged'] or info['residualMaxMm']>=.001:break
    x=z;final_rho=rho;multipliers=band.update(x)
    if row['outsideCount']==0:
        accepted=True
        break
if accepted:
    # Saved updated multipliers are the effective constraint forces at this solve.
    if multipliers is None:multipliers=np.zeros((2,len(tets)))
    np.savez_compressed(out/'volume-band-candidate.npz',nodes=x,rest=rest,tets=tets,faces=faces,fixed=fixed,classes=a['classes'],index=index,center=a['center'],multipliers=multipliers)
report={'index':index,'checkpoint':os.environ.get('QA_BAND_CHECKPOINT','bulk-3000.npz'),'finalRho':final_rho,'initial':initial,'linearizedFeasibility':diagnosis,'rows':rows,'localVolumeGuardPassed':accepted,'bulk':bulk,'targetBand':[.955,1.045],'engineeringGuard':[.95,1.05],'weightFloor':'1 / element count, only inequality term','productAccepted':False,'contactEvaluated':False,'selfIntersectionEvaluated':False,'physicalMaterialCalibrated':False}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print('LEVATOR_BAND_REPORT '+json.dumps(report),flush=True)
if not accepted:raise SystemExit(2)
