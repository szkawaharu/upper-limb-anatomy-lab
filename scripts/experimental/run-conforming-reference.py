"""Stateless source->target Dirichlet continuation, coarse material embedding.
Numerical loading path is not measured muscle activity or measured joint motion.
"""
import sys,json,time,hashlib,shutil
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from volumetric_reference import Solid
from cage_solver import minimize
from continuation import advance

p=Path(sys.argv[1]);elev=float(sys.argv[2]);plane=float(sys.argv[3]);bulk=float(sys.argv[4]if len(sys.argv)>4 else 30)
a=np.load(p/'cage.npz');refined=np.load(p/'conforming.npz');rest=a['nodes'];source=refined['points'];faces=refined['faces'];classes=refined['classes'];fixed=np.flatnonzero(classes>=0);solid=Solid(rest,a['tets'],fixed,element_weights=a['weights']);target=next(x for x in json.loads((p/'targets.json').read_text())['poses']if x['elev']==elev and x['plane']==plane);scale=float(a['worldUnitsPerMm']);center=a['centerWorld'];clock=time.time()
code_files=['run-conforming-reference.py','continuation.py','tetrahedral_cage.py','conforming_surface.py','cage_solver.py','projected_solid_metric.py','constitutive_metric.py','volumetric_reference.py'];code_hashes={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()for name in code_files};snapshot=p/'code-snapshot';snapshot.mkdir(exist_ok=True)
input_hashes={name:hashlib.sha256((p/name).read_bytes()).hexdigest()for name in ['cage.npz','conforming.npz','targets.json','preparation.json']}
for name in code_files:shutil.copy2(Path(__file__).with_name(name),snapshot/name)
RA=Rotation.from_quat(target['a']['q']).as_matrix();ta=np.array(target['a']['p'])-RA@np.array(target['a']['rest']);qa=np.array(target['a']['q']);qa[:3]*=-1;qb=np.array(target['b']['q']);q=np.r_[qa[3]*qb[:3]+qb[3]*qa[:3]+np.cross(qa[:3],qb[:3]),qa[3]*qb[3]-np.dot(qa[:3],qb[:3])];s=np.linalg.norm(q[:3]);angle=2*np.arctan2(s,q[3]);axis=q[:3]/s if s>1e-12 else np.zeros(3);anchor=(np.array(target['b']['rest'])-center)/scale;delta=(RA.T@(np.array(target['b']['p'])-ta)-np.array(target['b']['rest']))/scale

def transform(points,progress):return Rotation.from_rotvec(axis*angle*progress).apply(points-anchor)+anchor+progress*delta

def pins(progress):
 result=rest[fixed].copy();mask=classes[fixed]==1;result[mask]=transform(result[mask],progress);return result

def surface(nodes):return np.einsum('ni,nij->nj',refined['embeddingBarycentric'],nodes[refined['embeddingVertices']])

step=min(1.,np.deg2rad(3)/max(angle,1e-12));base_step=step;progress=0.;x=rest.copy();stages=[]
name='reference-e'+str(int(elev))+'-p'+str(int(plane))+'-k'+str(int(bulk));out=p/name
def solve(x,previous,next_progress):
 start=solid.predict_boundary(x,pins(next_progress))
 if not np.isfinite(solid.evaluate(start,bulk,gradient=False)[0]):
  return x,{'converged':False,'reason':'nonpositive-J-predictor','seconds':time.time()-clock}
 trial,report=minimize(solid,start,bulk=bulk,iterations=150)
 print(json.dumps({'fromProgress':previous,'attemptedProgress':next_progress,'seconds':time.time()-clock,**{k:v for k,v in report.items()if k!='rows'}}),flush=True)
 return trial,report
def save(x,progress,stages,rejections):
 now=surface(x);E=lambda points:np.stack([points[faces[:,k]]-points[faces[:,0]]for k in [1,2]],axis=2)
 R=E(source);N=E(now);area0=np.linalg.norm(np.cross(R[:,:,0],R[:,:,1]),axis=1);area1=np.linalg.norm(np.cross(N[:,:,0],N[:,:,1]),axis=1)
 edge=np.unique(np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0);er=np.linalg.norm(now[edge[:,1]]-now[edge[:,0]],axis=1)/np.linalg.norm(source[edge[:,1]]-source[edge[:,0]],axis=1)
 vol=lambda y:float(np.einsum('ti,ti->t',y[faces[:,0]],np.cross(y[faces[:,1]],y[faces[:,2]])).sum()/6)
 targets=source.copy();targets[refined['weights']==1]=transform(source[refined['weights']==1],progress);mask=np.isin(refined['weights'],[0,1]);error=float(np.linalg.norm(now[mask]-targets[mask],axis=1).max())
 F=(solid.edges(x)@solid.inverse)[refined['elements']];n0=np.cross(R[:,:,0],R[:,:,1]);expected=np.linalg.det(F)[:,None]*np.einsum('tij,tj->ti',np.linalg.inv(F).transpose(0,2,1),n0);actual=np.cross(N[:,:,0],N[:,:,1]);opposed=np.einsum('ti,ti->t',expected,actual)<0
 result={'sourceOpposedFaces':int(opposed.sum()),'experimental':True,'accepted':False,'target':target,'progress':progress,'bulk':bulk,'codeSha256':code_hashes,'cageNodes':len(rest),'sourceVertices':len(source),'originalSourceVertices':len(a['sourcePoints']),'surfaceVolumeRatio':vol(now)/vol(source),'surfaceEdgeRatio':[float(er.min()),float(er.max())],'surfaceAreaRatio':[float((area1/area0).min()),float((area1/area0).max())],'attachmentMaxMm':error,'allStagesConverged':all(r['converged']for r in stages),'seconds':time.time()-clock,'stages':stages,'rejectedTrials':rejections,'boneContactEvaluated':False,'selfCollisionEvaluated':False,'activation':None,'force':None}
 # Atomic individual files plus a content hash: a concurrent audit can reject
 # a mixed old/new pair rather than reading a half-written ZIP snapshot.
 temporary=out.with_suffix('.tmp.npz');np.savez_compressed(temporary,nodes=x,referenceNodes=rest,cageTets=a['tets'],sourceWorld=(now*scale+center)@RA.T+ta,referenceWorld=source*scale+center,sourceFaces=faces,sourceWeights=refined['weights'],embeddingVertices=refined['embeddingVertices'],embeddingBarycentric=refined['embeddingBarycentric'],triangleElements=refined['elements'],progress=progress)
 result['inputSha256']=input_hashes;result['snapshotSha256']=hashlib.sha256(temporary.read_bytes()).hexdigest();temporary.replace(out.with_suffix('.npz'));temporary=out.with_suffix('.tmp.json');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(out.with_suffix('.json'))
x,run=advance(rest,solve,base_step,on_accept=save)
save(x,run['progress'],run['accepted'],run['rejected'])
print('FINISHED',json.dumps({'progress':run['progress'],'completed':run['completed'],'seconds':time.time()-clock}),flush=True)
if not run['completed']:raise RuntimeError('Target not reached; only last converged state saved')
