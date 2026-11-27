"""Stateless source->target Dirichlet continuation, coarse material embedding.
Numerical loading path is not measured muscle activity or measured joint motion.
"""
import sys,json,time,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from volumetric_reference import Solid
from cage_solver import minimize

p=Path(sys.argv[1]);elev=float(sys.argv[2]);plane=float(sys.argv[3]);bulk=float(sys.argv[4]if len(sys.argv)>4 else 30)
a=np.load(p/'cage.npz');rest=a['nodes'];source=a['sourcePoints'];faces=a['sourceFaces'];classes=a['classes'];fixed=np.flatnonzero(classes>=0);solid=Solid(rest,a['tets'],fixed,element_weights=a['weights']);target=next(x for x in json.loads((p/'targets.json').read_text())['poses']if x['elev']==elev and x['plane']==plane);scale=float(a['worldUnitsPerMm']);center=a['centerWorld'];clock=time.time()
RA=Rotation.from_quat(target['a']['q']).as_matrix();ta=np.array(target['a']['p'])-RA@np.array(target['a']['rest']);qa=np.array(target['a']['q']);qa[:3]*=-1;qb=np.array(target['b']['q']);q=np.r_[qa[3]*qb[:3]+qb[3]*qa[:3]+np.cross(qa[:3],qb[:3]),qa[3]*qb[3]-np.dot(qa[:3],qb[:3])];s=np.linalg.norm(q[:3]);angle=2*np.arctan2(s,q[3]);axis=q[:3]/s if s>1e-12 else np.zeros(3);anchor=(np.array(target['b']['rest'])-center)/scale;delta=(RA.T@(np.array(target['b']['p'])-ta)-np.array(target['b']['rest']))/scale

def transform(points,progress):return Rotation.from_rotvec(axis*angle*progress).apply(points-anchor)+anchor+progress*delta

def pins(progress):
 result=rest[fixed].copy();mask=classes[fixed]==1;result[mask]=transform(result[mask],progress);return result

def surface(nodes):return np.einsum('ni,nij->nj',a['embeddingBarycentric'],nodes[a['embeddingVertices']])

step=min(1.,np.deg2rad(3)/max(angle,1e-12));base_step=step;progress=0.;x=rest.copy();stages=[]
name='reference-e'+str(int(elev))+'-p'+str(int(plane))+'-k'+str(int(bulk));out=p/name
while progress<1-1e-12:
 next_progress=min(1.,progress+step);start=solid.predict_boundary(x,pins(next_progress))
 if not np.isfinite(solid.evaluate(start,bulk,gradient=False)[0]):
  step*=.5
  if step<1e-7:raise RuntimeError('No positive-J continuation step')
  continue
 x,report=minimize(solid,start,bulk=bulk,iterations=100);progress=next_progress;stages.append({'progress':progress,'seconds':time.time()-clock,**report});print(json.dumps({k:v for k,v in stages[-1].items()if k!='rows'}),flush=True)
 now=surface(x);E=lambda points:np.stack([points[faces[:,k]]-points[faces[:,0]]for k in [1,2]],axis=2)
 R=E(source);N=E(now);area0=np.linalg.norm(np.cross(R[:,:,0],R[:,:,1]),axis=1);area1=np.linalg.norm(np.cross(N[:,:,0],N[:,:,1]),axis=1)
 edge=np.unique(np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0);er=np.linalg.norm(now[edge[:,1]]-now[edge[:,0]],axis=1)/np.linalg.norm(source[edge[:,1]]-source[edge[:,0]],axis=1)
 vol=lambda y:float(np.einsum('ti,ti->t',y[faces[:,0]],np.cross(y[faces[:,1]],y[faces[:,2]])).sum()/6)
 targets=source.copy();targets[a['sourceWeights']==1]=transform(source[a['sourceWeights']==1],progress);mask=np.isin(a['sourceWeights'],[0,1]);error=float(np.linalg.norm(now[mask]-targets[mask],axis=1).max())
 result={'experimental':True,'accepted':False,'target':target,'progress':progress,'bulk':bulk,'cageNodes':len(rest),'sourceVertices':len(source),'surfaceVolumeRatio':vol(now)/vol(source),'surfaceEdgeRatio':[float(er.min()),float(er.max())],'surfaceAreaRatio':[float((area1/area0).min()),float((area1/area0).max())],'attachmentMaxMm':error,'allStagesConverged':all(r['converged']for r in stages),'seconds':time.time()-clock,'stages':stages,'boneContactEvaluated':False,'selfCollisionEvaluated':False,'activation':None,'force':None}
 out.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');np.savez_compressed(out.with_suffix('.npz'),nodes=x,referenceNodes=rest,cageTets=a['tets'],sourceWorld=(now*scale+center)@RA.T+ta,referenceWorld=source*scale+center,sourceFaces=faces,sourceWeights=a['sourceWeights'],embeddingVertices=a['embeddingVertices'],embeddingBarycentric=a['embeddingBarycentric'])
 if not report['converged']:raise RuntimeError('Cage stage did not converge; partial result saved')
 step=min(base_step,step*1.25)
print('FINISHED',json.dumps({k:v for k,v in result.items()if k!='stages'}),flush=True)
