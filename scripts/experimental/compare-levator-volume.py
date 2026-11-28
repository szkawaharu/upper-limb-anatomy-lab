"""Compare local volume at a frozen, converged six-degree five-bone pose.
Do not change attachment points, source geometry, or equilibrium tolerance.
Passive material-energy sensitivity, not measured muscle properties."""
import json, os, time
from pathlib import Path
import numpy as np
from volumetric_reference import Solid
from cage_solver import minimize

src=Path(os.environ['QA_INPUT'])
out=Path(os.environ['QA_OUTPUT']);out.mkdir(parents=True,exist_ok=True)
files=list(src.glob('accepted-*-index-16.npz'))
assert len(files)==1,files
a=np.load(files[0]);assert int(a['index'])==16
rest=a['rest'];tets=a['tets'];fixed=a['fixed'];classes=a['classes'];x=a['nodes']
solid=Solid(rest,tets,fixed);pins=x[fixed].copy()
f=a['faces'];groups=classes[tets]
def audit(z,bulk):
    E,g,J=solid.evaluate(z,bulk)
    assert np.isfinite(E) and J.min()>1e-9
    log=np.log(J);absolute=np.abs(log)
    fully=np.all(groups>=0,axis=1)
    freecount=np.sum(groups<0,axis=1)
    def region(mask):
        return {'count':int(mask.sum()),'maxAbsLogJ':float(absolute[mask].max()) if mask.any() else None,
                'weightedMeanAbsLogJ':float(np.average(absolute[mask],weights=solid.volumes[mask])) if mask.any() else None}
    worst=[]
    for i in np.argsort(absolute)[-10:][::-1]:
        worst.append({'tetId':int(i),'J':float(J[i]),'referenceVolumeMm3':float(solid.volumes[i]),'fixedLabels':groups[i].tolist(),'condition':float(np.linalg.cond(solid.Dm[i])),'centroidAtlasMm':(rest[tets[i]].mean(axis=0)+a['center']).tolist()})
    p=rest[f];q=z[f]
    area=lambda v:np.linalg.norm(np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0]),axis=1)
    ratios=area(q)/area(p)
    return {'bulk':bulk,'JRange':[float(J.min()),float(J.max())],'maxAbsLogJ':float(absolute.max()),
      'volumeRatio':float(solid.weights@J),'volumeWeightedRmsLogJ':float(np.sqrt(solid.weights@(log*log))),
      'volumeFractionOutside5Percent':float(solid.weights@((J<.95)|(J>1.05))),
      'byFreeVertices':{str(i):region(freecount==i) for i in range(5)},
      'fullyFixedMixedBoneElements':int(np.sum(fully & np.any(groups!=groups[:,:1],axis=1))),
      'surfaceAreaRatioRange':[float(ratios.min()),float(ratios.max())],
      'residualMm':float(np.linalg.norm(solid.factor.solve(g[solid.free]),axis=1).max()),
      'attachmentMaxMm':float(np.linalg.norm(z[fixed]-pins,axis=1).max()),'worstElements':worst}
rows=[{'converged':True,**audit(x,30)}]
assert rows[0]['residualMm']<.001
for bulk in [100,300,1000]:
    start=time.monotonic()
    z,info=minimize(solid,x,bulk=bulk,iterations=30,tolerance_mm=.001)
    row={**audit(z,bulk),'converged':info['converged'],'reason':info['reason'],'iterations':info['iterations'],'seconds':time.monotonic()-start}
    assert row['attachmentMaxMm']==0
    rows.append(row)
    print('LEVATOR_BULK_CASE '+json.dumps(row),flush=True)
    if info['converged']:
        np.savez_compressed(out/f'bulk-{bulk}.npz',nodes=z,rest=rest,tets=tets,faces=f,fixed=fixed,classes=classes,index=16,center=a['center'])
        x=z
    else:
        break
report={'sourceCheckpoint':files[0].name,'rows':rows,'productAccepted':False,'anatomicalAcceptance':False,'contactEvaluated':False,'selfIntersectionEvaluated':False,'physicalMaterialCalibrated':False}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print('LEVATOR_BULK_REPORT '+json.dumps(report),flush=True)
if not all(row['converged'] for row in rows):raise SystemExit(2)
