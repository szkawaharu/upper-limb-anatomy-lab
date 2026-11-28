"""Actual-frame continuation with local inequality constraints at each pose.
The last accepted constrained shape seeds the next predictor. Multipliers are
local to each trial; no rejected state or multiplier leaks into the path.
Offline numerical candidate only, not a calibrated material or release gate.
"""
import os,json,time,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from volumetric_reference import Solid
from volumetric_preconditioner import coupled_boundary_predictor
from local_volume_band import VolumeBand
from local_volume_constraints import log_volume_jacobian
from contact_potential import AugmentedSolid
from cage_solver import minimize
from joint_path import follow

src=Path(os.environ['QA_INPUT']);out=Path(os.environ['QA_OUTPUT']);out.mkdir(parents=True,exist_ok=True)
files=list(src.glob('accepted-*-index-16.npz'));assert files
a=np.load(sorted(files)[0])
rest=a['rest'];tets=a['tets'];fixed=a['fixed'];classes=a['classes'];faces=a['faces'];center=a['center']
assert np.array_equal(fixed,np.flatnonzero(classes>=0))
assert np.all(np.isin(fixed,np.unique(faces))) and set(classes[fixed])==set(range(5))
solid=Solid(rest,tets,fixed)
bulk=1000.
framefile=Path(os.environ.get('QA_FRAME_FILE',str(src/'levator-solid-frames.json')))
frames=json.loads(framefile.read_text())['frames']
maximum=int(os.environ.get('QA_MAX_INDEX','40'))
assert maximum in (40,80,160) and len(frames)>maximum
assert all(f['index']==i and f['elev']==i*.375 and f['neck']==i*.375 and f['plane']==40 and f['backend']==('qa-additive-preview' if f['elev']>40 else 'validated-field') for i,f in enumerate(frames))
# Verify extended grid against the exact original frames on the shared region.
old=json.loads((src/'levator-solid-frames.json').read_text())['frames']
for i,f in enumerate(old):
    for key in ['bones','sourceOrigin','scale']:
        if key=='bones':
            for bone in f[key]:
                for k in ['q','p','rest']:
                    assert np.allclose(f[key][bone][k],frames[i][key][bone][k],rtol=0,atol=1e-12)
        else:assert np.allclose(f[key],frames[i][key],rtol=0,atol=1e-12)
Q=np.array([[-1,0,0],[0,0,1],[0,1,0]],float);keys=['C1','C2','C3','C4','scapula']
def prescribed(frame):
    o=np.asarray(frame['sourceOrigin']);s=frame['scale'];shift=np.array([0,1.70,.08])
    world=(rest[fixed]+center-o)@Q.T*s+shift
    moved=np.empty_like(world)
    for k,key in enumerate(keys):
        mask=classes[fixed]==k;b=frame['bones'][key]
        moved[mask]=Rotation.from_quat(b['q']).apply(world[mask]-b['rest'])+b['p']
    return ((moved-shift)/s)@Q+o-center
targets=[prescribed(f) for f in frames]
assert np.max(np.linalg.norm(targets[0]-rest[fixed],axis=1))<1e-7
tri=rest[faces];area0=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
edges=np.unique(np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0)
edge0=np.linalg.norm(rest[edges[:,1]]-rest[edges[:,0]],axis=1)
assert area0.min()>0 and edge0.min()>0
def metrics(z):
    E,g,J=solid.evaluate(z,bulk);assert np.isfinite(E) and np.isfinite(J).all() and J.min()>1e-9
    tri=z[faces];ar=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/area0
    er=np.linalg.norm(z[edges[:,1]]-z[edges[:,0]],axis=1)/edge0
    return {'JRange':[float(J.min()),float(J.max())],'outsideCount':int(((J<.95)|(J>1.05)).sum()),'volumeRatio':float(solid.weights@J),'surfaceAreaRatioRange':[float(ar.min()),float(ar.max())],'edgeRatioRange':[float(er.min()),float(er.max())],'baseResidualMm':float(np.linalg.norm(solid.factor.solve(g[solid.free]),axis=1).max()),'nonPositiveElements':int((J<=0).sum())}
destinations=list(range(8,maximum+1,8))+list(range(maximum-8,-1,-8))
start_index=0;start_nodes=rest.copy()
effective=None
seed_file=os.environ.get('QA_CONSTRAINED_SEED')
if seed_file:
    seed=np.load(seed_file)
    for k,v in [('rest',rest),('tets',tets),('fixed',fixed),('classes',classes),('faces',faces),('center',center)]:
        assert np.array_equal(seed[k],v),k
    start_index=int(seed['index']);assert start_index==80 and maximum==160
    start_nodes=seed['nodes'].copy();effective=seed['multipliers'].copy()
    assert effective.shape==(2,len(tets)) and np.isfinite(effective).all() and np.all(effective>=0)
    assert np.max(np.linalg.norm(start_nodes[fixed]-targets[start_index],axis=1))<1e-9
    E,g,J=solid.evaluate(start_nodes,bulk)
    assert np.isfinite(E) and J.min()>=.95 and J.max()<=1.05
    band=VolumeBand(solid,30)
    _,A=log_volume_jacobian(solid,start_nodes,band.cells)
    g+=(A.T@(band.weights*(effective[1]-effective[0]))).reshape(start_nodes.shape)
    residual=float(np.linalg.norm(solid.factor.solve(g[solid.free]),axis=1).max())
    assert residual<.001
    destinations=list(range(start_index+8,maximum+1,8))+list(range(maximum-8,-1,-8))
    np.savez_compressed(out/'accepted-000-index-80.npz',nodes=start_nodes,rest=rest,tets=tets,faces=faces,fixed=fixed,classes=classes,index=start_index,center=center,multipliers=effective)
    print('LEVATOR_CONSTRAINED_SEED '+json.dumps({'index':start_index,'residualMm':residual,'JRange':[float(J.min()),float(J.max())]}),flush=True)
def solve_at(x,previous,index):
    global effective
    start=time.monotonic();target=targets[index]
    z,predictor=coupled_boundary_predictor(solid,x,target,bulk)
    _,_,J=solid.evaluate(z,bulk,False)
    base={'elev':frames[index]['elev'],'neck':frames[index]['neck'],'predictor':predictor,'predictorJRange':[float(J.min()),float(J.max())]}
    if not np.isfinite(J).all() or J.min()<=1e-9:
        return x,{**base,'converged':False,'reason':'nonpositive-predictor'}
    rows=[];multipliers=None
    for rho in [30,100,300,1000,3000,10000]:
        band=VolumeBand(solid,rho,multipliers=multipliers)
        augmented=AugmentedSolid(solid,[band])
        trial,info=minimize(augmented,z,bulk=bulk,iterations=40,tolerance_mm=.001)
        m=metrics(trial);pin=float(np.linalg.norm(trial[fixed]-target,axis=1).max())
        # Recompute equilibrium using effective multipliers at the final point.
        _,g,_=solid.evaluate(trial,bulk)
        _,A=log_volume_jacobian(solid,trial,band.cells)
        p=band.update(trial)
        g+=(A.T@(band.weights*(p[1]-p[0]))).reshape(trial.shape)
        residual=float(np.linalg.norm(solid.factor.solve(g[solid.free]),axis=1).max())
        row={'rho':rho,'iterations':info['iterations'],'innerConverged':info['converged'],'augmentedResidualMm':residual,'attachmentErrorMm':pin,**m}
        rows.append(row)
        if not info['converged'] or residual>=.001 or pin>=1e-9:break
        z=trial;multipliers=p
        if m['outsideCount']==0:
            effective=p.copy()
            result={**base,**row,'converged':True,'reason':'equilibrium-and-local-volume-guard','stages':rows,'maxVertexStepMm':float(np.linalg.norm(z-x,axis=1).max()),'seconds':time.monotonic()-start}
            print('LEVATOR_CONSTRAINED_TRIAL '+json.dumps(result),flush=True)
            return z,result
    result={**base,'converged':False,'reason':'local-band-or-equilibrium-failed','stages':rows,'seconds':time.monotonic()-start}
    print('LEVATOR_CONSTRAINED_TRIAL '+json.dumps(result),flush=True)
    return x,result
def accepted(x,index,records,rejected):
    assert effective is not None
    np.savez_compressed(out/f'accepted-{len(records):03d}-index-{index:02d}.npz',nodes=x,rest=rest,tets=tets,faces=faces,fixed=fixed,classes=classes,index=index,center=center,multipliers=effective)
    (out/'path-progress.json').write_text(json.dumps({'index':index,'accepted':records,'rejected':rejected},indent=2)+'\n')
(out/'levator-solid-frames.json').write_bytes(framefile.read_bytes())
x,path=follow(start_nodes,start_index,destinations,solve_at,accepted)
report={'requestedElevationsDeg':[start_index*.375,*[i*.375 for i in destinations]],'seedFile':seed_file,'frameSHA256':hashlib.sha256(framefile.read_bytes()).hexdigest(),'bulk':bulk,'path':path,'allAcceptedLocalGuardPassed':bool(path['accepted']) and all(v['outsideCount']==0 for v in path['accepted']),'neutralReturnMaxMm':float(np.linalg.norm(x-rest,axis=1).max()) if path['completed'] else None,'productAccepted':False,'contactEvaluated':False,'continuousCollisionCertified':False,'physicalMaterialCalibrated':False}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print('LEVATOR_CONSTRAINED_PATH '+json.dumps(report),flush=True)
if not path['completed']:raise SystemExit(2)
