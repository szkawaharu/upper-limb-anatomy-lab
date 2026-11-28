"""Fixed-pose untangling experiment on audited stretch candidates.
Keep original five-bone pins and local volume guard. A zero result from the
independent finite-edge crossing screen is necessary, not sufficient for
anatomical or continuous-contact acceptance."""
import os,json,time
from pathlib import Path
import numpy as np
from volumetric_reference import Solid
from local_volume_band import VolumeBand
from local_stretch_band import StretchBand
from local_volume_constraints import log_volume_jacobian
from local_self_separation import SelfSeparation
from surface_crossing_screen import screen
from contact_potential import AugmentedSolid
from cage_solver import minimize

src=Path(os.environ['QA_INPUT']);out=Path(os.environ['QA_OUTPUT']);out.mkdir(parents=True,exist_ok=True)
a=np.load(src/'volume-band-candidate.npz')
reference=np.load(Path(os.environ['QA_REFERENCE']))
assert int(reference['index'])==80
rest=a['rest'];tets=a['tets'];faces=a['faces'];fixed=a['fixed']
for k in ['rest','tets','faces','fixed','classes','center']:
    assert np.array_equal(a[k],reference[k]),k
ref=reference['nodes'];start_nodes=a['nodes'].copy();x=start_nodes.copy();pins=x[fixed].copy()
solid=Solid(rest,tets,fixed);bulk=1000.;strength=float(os.environ['QA_STRETCH_STRENGTH'])
stretch=StretchBand(solid,strength)
multipliers=a['multipliers'].copy()
assert multipliers.shape==(2,len(tets)) and np.all(np.isfinite(multipliers)) and np.all(multipliers>=0)
ref_cross,_=screen(ref,faces);assert not ref_cross
initial_cross,_=screen(x,faces)
assert initial_cross
edges=np.unique(np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0)
def areas(v):
    p=v[faces];return np.linalg.norm(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]),axis=1)
area0=areas(rest);assert area0.min()>0
edge0=np.linalg.norm(rest[edges[:,1]]-rest[edges[:,0]],axis=1)
def metrics(v):
    E,g,J=solid.evaluate(v,bulk);assert np.isfinite(E) and J.min()>1e-9
    s=np.linalg.svd(solid.edges(v)@solid.inverse,compute_uv=False)
    ar=areas(v)/area0;er=np.linalg.norm(v[edges[:,1]]-v[edges[:,0]],axis=1)/edge0
    return {'JRange':[float(J.min()),float(J.max())],'outsideVolumeCount':int(((J<.95)|(J>1.05)).sum()),'principalStretchRange':[float(s.min()),float(s.max())],'maxConditionF':float(np.max(s[:,0]/s[:,-1])),'areaRatioRange':[float(ar.min()),float(ar.max())],'edgeRatioRange':[float(er.min()),float(er.max())],'attachmentErrorMm':float(np.linalg.norm(v[fixed]-pins,axis=1).max())}
initial=metrics(x);assert initial['outsideVolumeCount']==0
# Recheck the saved constrained equilibrium with its original stretch term.
band=VolumeBand(solid,300)
_,A=log_volume_jacobian(solid,x,band.cells);_,g,_=solid.evaluate(x,bulk)
g+=(A.T@(band.weights*(multipliers[1]-multipliers[0]))).reshape(x.shape)
g+=stretch.evaluate(x)[1]
seed_res=float(np.linalg.norm(solid.factor.solve(g[solid.free]),axis=1).max());assert seed_res<.001
pairs=set(initial_cross);rows=[];accepted=False;reason='schedule-exhausted'
for rho_contact in [.01,.1,1.,10.,100.]:
    # Keep the reference-side orientation fixed relative to the original
    # target pose. Add newly detected pairs, not ad hoc vertex displacements.
    try:contact=SelfSeparation(solid,faces,ref,start_nodes,pairs,rho_contact)
    except ValueError as exc:
        reason=str(exc);break
    if contact.fixedConflictCount:
        reason='fixed-separation-conflict';break
    z=x.copy();p=multipliers.copy();stages=[];candidate=None
    for rho_volume in [300,1000,3000,10000]:
        start=time.monotonic();band=VolumeBand(solid,rho_volume,multipliers=p)
        aug=AugmentedSolid(solid,[band,stretch,contact])
        trial,info=minimize(aug,z,bulk=bulk,iterations=40,tolerance_mm=.001)
        assert np.array_equal(trial[fixed],pins)
        _,g,_=aug.evaluate(trial,bulk)
        residual=float(np.linalg.norm(solid.factor.solve(g[solid.free]),axis=1).max())
        m=metrics(trial);violation=contact.violation(trial)
        row={'contactRho':rho_contact,'volumeRho':rho_volume,'converged':bool(info['converged'] and residual<.001),'residualMm':residual,'iterations':info['iterations'],'pairCount':len(pairs),'maxSeparationViolationMm':float(violation.max(initial=0)),'seconds':time.monotonic()-start,**m}
        stages.append(row);print('LEVATOR_SELF_SEPARATION_STAGE '+json.dumps(row),flush=True)
        if not row['converged']:break
        z=trial;p=band.update(z)
        if m['outsideVolumeCount']==0:
            candidate=z.copy();break
    if candidate is None:
        rows.append({'contactRho':rho_contact,'candidateSaved':False,'stages':stages});reason='volume-or-equilibrium-failed';break
    crossings,screen_info=screen(candidate,faces)
    row={'contactRho':rho_contact,'candidateSaved':True,'crossingCount':len(crossings),'crossingPairs':sorted(crossings),'screen':screen_info,'separation':contact.details,'stages':stages}
    rows.append(row)
    folder=out/f'contact-{rho_contact:g}';folder.mkdir(exist_ok=True)
    np.savez_compressed(folder/'volume-band-candidate.npz',nodes=candidate,rest=rest,tets=tets,faces=faces,fixed=fixed,classes=a['classes'],index=a['index'],center=a['center'],multipliers=p)
    print('LEVATOR_SELF_SEPARATION_CASE '+json.dumps({'contactRho':rho_contact,'crossingCount':len(crossings),'maxSeparationViolationMm':stages[-1]['maxSeparationViolationMm'],'JRange':stages[-1]['JRange'],'residualMm':stages[-1]['residualMm']}),flush=True)
    # Preserve every accepted numerical candidate, even if crossing screen fails.
    x=candidate;multipliers=p
    if not crossings and stages[-1]['maxSeparationViolationMm']<=1e-6:
        accepted=True;reason='sampled-separation-and-volume-guard';break
    pairs.update(crossings)
report={'index':int(a['index']),'stretchStrength':strength,'referenceIndex':80,'initial':initial,'seedResidualMm':seed_res,'initialCrossingPairs':sorted(initial_cross),'rows':rows,'sampledSeparationGuardPassed':accepted,'reason':reason,'productAccepted':False,'boneContactEvaluated':False,'continuousCollisionCertified':False,'coplanarOverlapCertified':False,'sharedVertexPairsExcluded':True,'physicalMaterialCalibrated':False}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print('LEVATOR_SELF_SEPARATION_REPORT '+json.dumps({**report,'rows':[{'contactRho':r['contactRho'],'candidateSaved':r['candidateSaved'],'crossingCount':r.get('crossingCount'),'stages':r['stages']} for r in rows]}),flush=True)
if not accepted:raise SystemExit(2)
