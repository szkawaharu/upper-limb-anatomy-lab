"""Offline five-bone material experiment on the exact-source principal shell.
The atlas-mm coordinate system is converted through actual browser frames.
Near-surface labels are experimental, not measured anatomical footprints.
No active force, contact, physical material constants, or product acceptance.
"""
import json,os,hashlib,urllib.request,sys,time
from pathlib import Path
import numpy as np
import trimesh
from scipy.spatial.transform import Rotation
from volumetric_reference import Solid
from joint_path import follow
from volumetric_gauss_newton import minimize as minimize_gn
from cage_solver import minimize as minimize_tangent
from volumetric_preconditioner import factor_rest_metric
from scipy.sparse.linalg import LinearOperator

root=Path('test-results/v076-levator-solid')
out=Path(os.environ.get('QA_OUTPUT','test-results/v076-levator-solid-path'));out.mkdir(parents=True,exist_ok=True)
tetra=Path(os.environ.get('QA_LEVATOR_TETRA_DIR',str(root/'tetra')))
meta=json.loads((tetra/'reference-mesh-audit.json').read_text())
assert meta['sourceBoundaryExact'] and meta['sourcePrefixMaxError']==0
shell=np.load(root/'main-shell.npz');mesh=np.load(tetra/'v58Levator-tetrahedra.npz')
nodes=np.asarray(mesh['nodes']);tets=mesh['tets'];faces=shell['faces']
assert np.array_equal(nodes[:len(shell['points'])],shell['points'])
assert meta['sourceSha256']==json.loads(Path('docs/validation/v076/levator-exact-source-solid.json').read_text())['tetrahedralization']['sourceSha256']
evidence=json.loads(Path('docs/validation/v076/levator-attachment-surfaces.json').read_text())
cache=Path(os.environ.get('ATLAS_CACHE','/tmp/upper-limb-atlas'));cache.mkdir(parents=True,exist_ok=True)
url='https://raw.githubusercontent.com/Kevin-Mattheus-Moerman/BodyParts3D/2c83c698788cff57f6317bddcf712e7acfbea7d8/assets/BodyParts3D_data/stl/'
files=['FMA12519.stl','FMA12520.stl','FMA12521.stl','FMA12522.stl','FMA13395.stl']
keys=['C1','C2','C3','C4','scapula']
p=shell['points'];dist=[]
for file in files:
    src=next(x for x in evidence['sources'] if x['file']==file);path=cache/file
    if not path.exists():
        with urllib.request.urlopen(url+file,timeout=90) as f:path.write_bytes(f.read())
    assert hashlib.sha256(path.read_bytes()).hexdigest()==src['sha256']
    bone=trimesh.load_mesh(path,process=True)
    dist.append(np.concatenate([trimesh.proximity.closest_point(bone,p[i:i+128])[1] for i in range(0,len(p),128)]))
d=np.stack(dist,axis=1);nearest=d.argmin(axis=1);gap=d.min(axis=1)
classes=np.full(len(nodes),-1,int);classes[:len(p)]=np.where(gap<=1.,nearest,-1)
fixed=np.flatnonzero(classes>=0);counts=np.bincount(classes[fixed],minlength=5)
assert np.all(counts>0)
# All declared attachment nodes must be source boundary vertices.
assert np.all(np.isin(fixed,np.unique(faces)))
np.savez_compressed(out/'attachments.npz',classes=classes,sourceDistancesMm=d,fixed=fixed)
report={'radiusMm':1,'counts':dict(zip(keys,counts.tolist())),'fixed':len(fixed),'nodes':len(nodes),'tets':len(tets),'interiorNodes':int(len(nodes)-len(np.unique(faces))),'sourceSha256':meta['sourceSha256'],'productAccepted':False}
print('LEVATOR_SOLID_ATTACHMENTS '+json.dumps(report),flush=True)
frames=json.loads((out/'levator-solid-frames.json').read_text())['frames']
assert 17<=len(frames)<=107 and all(x['index']==i and x['elev']==i*.375 and x['neck']==i*.375 and x['plane']==40 and x['backend']=='validated-field' for i,x in enumerate(frames))
# Work in centered atlas millimetres, not browser units. World chart x=-atlasX,
# y=atlasZ, z=atlasY; origin offsets are accounted for explicitly.
center=nodes.mean(axis=0);rest=nodes-center
Q=np.array([[-1,0,0],[0,0,1],[0,1,0]],float)
def prescribed(frame):
    s=frame['scale'];o=np.asarray(frame['sourceOrigin']);shift=np.array([0,1.70,.08])
    world=(nodes[fixed]-o)@Q.T*s+shift
    moved=np.empty_like(world)
    for k,key in enumerate(keys):
        mask=classes[fixed]==k;b=frame['bones'][key]
        moved[mask]=Rotation.from_quat(b['q']).apply(world[mask]-b['rest'])+b['p']
    return ((moved-shift)/s)@Q+o-center
targets=[prescribed(f) for f in frames]
assert np.max(np.linalg.norm(targets[0]-rest[fixed],axis=1))<1e-7,'Nonidentity neutral boundary'
solid=Solid(rest,tets,fixed)
bulk=float(os.environ.get('QA_LEVATOR_BULK','30'))
iterations=int(os.environ.get('QA_LEVATOR_ITERATIONS','150'))
solver=os.environ.get('QA_LEVATOR_SOLVER','lbfgs')
assert solver in ('lbfgs','coupled-gn','exact-tangent')
preconditioner=None;preconditioner_info=None
if solver=='coupled-gn':
    start=time.monotonic();factor,preconditioner_info=factor_rest_metric(solid,bulk)
    n=len(solid.free)*3
    preconditioner=LinearOperator((n,n),matvec=factor.solve,dtype=float)
    preconditioner_info['setupSeconds']=time.monotonic()-start
    print('LEVATOR_SOLID_PRECONDITIONER '+json.dumps(preconditioner_info),flush=True)
_,_,J=solid.evaluate(rest,bulk);assert np.max(abs(J-1))<1e-9
start_index=0;start_nodes=rest.copy()
destinations=[8,16,8,0]
seed_file=os.environ.get('QA_LEVATOR_SEED')
if seed_file:
    seed=np.load(seed_file)
    assert np.array_equal(seed['rest'],rest) and np.array_equal(seed['tets'],tets) and np.array_equal(seed['faces'],faces)
    assert np.array_equal(seed['fixed'],fixed) and np.array_equal(seed['classes'],classes) and np.array_equal(seed['center'],center)
    start_index=int(seed['index']);assert start_index==16
    start_nodes=seed['nodes'].copy()
    assert np.max(np.linalg.norm(start_nodes[fixed]-targets[start_index],axis=1))<1e-9
    E,g,J=solid.evaluate(start_nodes,bulk)
    residual=float(np.linalg.norm(solid.factor.solve(g[solid.free]),axis=1).max())
    assert np.isfinite(E) and J.min()>1e-9 and residual<.001
    assert J.min()>=.95 and J.max()<=1.05,'Seed fails experimental local volume guard'
    destinations=[24,32,40,32,24,16,8,0]
    assert len(frames)>max(destinations)
    np.savez_compressed(out/'accepted-000-index-16.npz',nodes=start_nodes,rest=rest,tets=tets,faces=faces,fixed=fixed,classes=classes,index=start_index,center=center)
    print('LEVATOR_EXTENDED_SEED '+json.dumps({'index':start_index,'residualMm':residual,'JRange':[float(J.min()),float(J.max())],'bulk':bulk}),flush=True)
rows=[]
def solve_at(x,previous,index):
    t=time.monotonic();target=targets[index]
    initial=solid.predict_boundary(x,target)
    _,_,predJ=solid.evaluate(initial,bulk,False)
    base={'elev':frames[index]['elev'],'neck':frames[index]['neck'],'predictorMinJ':float(predJ.min()),'predictorMaxJ':float(predJ.max()),'fixedStepMm':float(np.linalg.norm(target-targets[previous],axis=1).max())}
    if not np.isfinite(predJ).all() or predJ.min()<=1e-9:
        result={**base,'converged':False,'reason':'nonpositive-predictor','seconds':time.monotonic()-t}
        print('LEVATOR_SOLID_TRIAL '+json.dumps(result),flush=True);return x,result
    if solver=='coupled-gn':
        z,info=minimize_gn(solid,initial,bulk=bulk,iterations=iterations,residual_tolerance_mm=.001,preconditioner=preconditioner,preconditioner_info=preconditioner_info,cg_iterations=150,cg_rtol=1e-3)
        info['linearIterationsTotal']=sum(row['cgIterations'] for row in info.pop('rows'))
    elif solver=='exact-tangent':
        z,info=minimize_tangent(solid,initial,bulk=bulk,iterations=iterations,tolerance_mm=.001)
        info['tangentSteps']=info.pop('rows')
        info['preconditionedGradientMaxMm']=info['residualMaxMm']
        info['residualToleranceMm']=.001
    else:
        z,info=solid.minimize(initial,bulk=bulk,iterations=iterations,residual_tolerance_mm=.001)
    _,_,J=solid.evaluate(z,bulk,False)
    if not np.isfinite(J).all() or np.min(J)<=1e-9:raise RuntimeError('Solver returned invalid material state')
    extrema=[]
    for label,k in [('minJ',int(J.argmin())),('maxJ',int(J.argmax()))]:
        ids=tets[k];extrema.append({'kind':label,'tetId':k,'J':float(J[k]),'sourceCentroidAtlasMm':nodes[ids].mean(axis=0).tolist(),'referenceVolumeMm3':float(solid.volumes[k]),'referenceCondition':float(np.linalg.cond(solid.Dm[k])),'fixedLabels':classes[ids].tolist(),'vertexIds':ids.tolist()})
    info={**base,**info,'minJ':float(J.min()),'maxJ':float(J.max()),'JQuantiles':np.quantile(J,[0,.01,.5,.99,1]).tolist(),'volumeRatio':float(np.dot(solid.weights,J)),'nonPositiveElements':int((J<=0).sum()),'extrema':extrema,'fullyFixedElements':int(np.all(classes[tets]>=0,axis=1).sum()),'attachmentMaxMm':float(np.linalg.norm(z[fixed]-target,axis=1).max()),'maxAbsLogJ':float(np.max(np.abs(np.log(J)))),'maxVertexStepMm':float(np.linalg.norm(z-x,axis=1).max()),'seconds':time.monotonic()-t}
    assert info['attachmentMaxMm']<1e-9
    print('LEVATOR_SOLID_TRIAL '+json.dumps(info),flush=True)
    return z,info
def accepted(x,index,records,rejected):
    np.savez_compressed(out/'last-accepted.npz',nodes=x,rest=rest,tets=tets,faces=faces,fixed=fixed,classes=classes,index=index,center=center)
    if index in (0,8,16,24,32,40):
        np.savez_compressed(out/f'accepted-{len(records):03d}-index-{index:02d}.npz',nodes=x,rest=rest,tets=tets,faces=faces,fixed=fixed,classes=classes,index=index,center=center)
    (out/'path-progress.json').write_text(json.dumps({'index':index,'accepted':records,'rejected':rejected},indent=2)+'\n')
x,path=follow(start_nodes,start_index,destinations,solve_at,accepted)
summary={**report,'solver':solver,'tetraInputSha256':meta['inputSha256'],'bulkPenalty':bulk,'requestedElevationsDeg':[i*.375 for i in [start_index,*destinations]],'requestedNeckDeg':[i*.375 for i in [start_index,*destinations]],'seedFile':seed_file,'experimentalLocalVolumeGuard':[.95,1.05],'allAcceptedWithinLocalVolumeGuard':bool(path['accepted']) and all(v['minJ']>=.95 and v['maxJ']<=1.05 for v in path['accepted']),'path':path,'neutralReturnMaxMm':float(np.linalg.norm(x-rest,axis=1).max()) if path['completed'] else None,'anatomicalAcceptance':False,'physicalVolumeVerified':False,'contactEvaluated':False,'limitations':['principal source shell only','unrefined tetrahedra; all source nodes on boundary' if report['interiorNodes']==0 else 'interior-refined tetrahedra; source boundary unchanged','near-surface masks not measured anatomical footprints','normalized passive energy; no active force or bone contact']}
(out/'report.json').write_text(json.dumps(summary,indent=2)+'\n')
print('LEVATOR_SOLID_PATH '+json.dumps(summary),flush=True)
# A nonconverged attempt remains a measured failure, never a completed path.
if not path['completed']:sys.exit(2)
