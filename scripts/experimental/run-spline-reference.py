"""Stateless source->target Dirichlet continuation, coarse material embedding.
Numerical loading path is not measured muscle activity or measured joint motion.
"""
import sys,json,time,hashlib,shutil
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from spline_cage import SplineCage,SplineSolid

p=Path(sys.argv[1]);elev=float(sys.argv[2]);plane=float(sys.argv[3]);bulk=float(sys.argv[4]if len(sys.argv)>4 else 30)
iterations=int(sys.argv[5]if len(sys.argv)>5 else 200)
a=np.load(p/'cage.npz');rest=a['nodes'];source=a['sourcePoints'];faces=a['sourceFaces'];classes=a['classes'];fixed=np.flatnonzero(classes>=0);c=SplineCage(source,float(a['spacingMm']));assert np.array_equal(rest,c.nodes);solid=SplineSolid(c,a['quadraturePoints'],a['quadratureWeights'],classes);source_embedding=c.embed(source);solid.guards=c.embed(np.concatenate([source,source[faces].mean(1)]));target=next(x for x in json.loads((p/'targets.json').read_text())['poses']if x['elev']==elev and x['plane']==plane);scale=float(a['worldUnitsPerMm']);center=a['centerWorld'];clock=time.time()
code_files=['run-spline-reference.py','spline_cage.py','spline_regularization.py','constitutive_metric.py'];code_hashes={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()for name in code_files};snapshot=p/'code-snapshot';snapshot.mkdir(exist_ok=True)
for name in code_files:shutil.copy2(Path(__file__).with_name(name),snapshot/name)
RA=Rotation.from_quat(target['a']['q']).as_matrix();ta=np.array(target['a']['p'])-RA@np.array(target['a']['rest']);qa=np.array(target['a']['q']);qa[:3]*=-1;qb=np.array(target['b']['q']);q=np.r_[qa[3]*qb[:3]+qb[3]*qa[:3]+np.cross(qa[:3],qb[:3]),qa[3]*qb[3]-np.dot(qa[:3],qb[:3])];s=np.linalg.norm(q[:3]);angle=2*np.arctan2(s,q[3]);axis=q[:3]/s if s>1e-12 else np.zeros(3);anchor=(np.array(target['b']['rest'])-center)/scale;delta=(RA.T@(np.array(target['b']['p'])-ta)-np.array(target['b']['rest']))/scale

def transform(points,progress):return Rotation.from_rotvec(axis*angle*progress).apply(points-anchor)+anchor+progress*delta

def pins(progress):
 result=rest[fixed].copy();mask=classes[fixed]==1;result[mask]=transform(result[mask],progress);return result

def surface(nodes):return c.positions(nodes,source_embedding)

step=min(1.,np.deg2rad(3)/max(angle,1e-12));base_step=step;progress=0.;x=rest.copy();stages=[]
name='reference-e'+str(int(elev))+'-p'+str(int(plane))+'-k'+str(int(bulk));out=p/name
while progress<1-1e-12:
 next_progress=min(1.,progress+step);start=solid.predict(x,pins(next_progress))
 if not np.isfinite(solid.evaluate(start,bulk)[0]):
  step*=.5
  if step<1e-7:raise RuntimeError('No positive-J continuation step')
  continue
 x,report=solid.minimize(start,bulk=bulk,iterations=iterations,callback=lambda row:print(json.dumps({'event':'iteration','progress':next_progress,'seconds':time.time()-clock,**row}),flush=True)if row['iteration']==1 or row['iteration']%10==0 else None);progress=next_progress;stages.append({'progress':progress,'seconds':time.time()-clock,**report});print(json.dumps({k:v for k,v in stages[-1].items()if k!='rows'}),flush=True)
 now=surface(x);E=lambda points:np.stack([points[faces[:,k]]-points[faces[:,0]]for k in [1,2]],axis=2)
 R=E(source);N=E(now);area0=np.linalg.norm(np.cross(R[:,:,0],R[:,:,1]),axis=1);area1=np.linalg.norm(np.cross(N[:,:,0],N[:,:,1]),axis=1)
 edge=np.unique(np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0);er=np.linalg.norm(now[edge[:,1]]-now[edge[:,0]],axis=1)/np.linalg.norm(source[edge[:,1]]-source[edge[:,0]],axis=1)
 vol=lambda y:float(np.einsum('ti,ti->t',y[faces[:,0]],np.cross(y[faces[:,1]],y[faces[:,2]])).sum()/6)
 targets=source.copy();targets[a['sourceWeights']==1]=transform(source[a['sourceWeights']==1],progress);mask=np.isin(a['sourceWeights'],[0,1]);error=float(np.linalg.norm(now[mask]-targets[mask],axis=1).max())
 guardJ=np.linalg.det(c.gradients(x,solid.guards));F=c.gradients(x,c.embed(source[faces].mean(1)));normal=np.cross(R[:,:,0],R[:,:,1]);mapped=np.linalg.det(F)[:,None]*np.einsum('tij,tj->ti',np.linalg.inv(F).transpose(0,2,1),normal);opposed=np.einsum('ti,ti->t',mapped,np.cross(N[:,:,0],N[:,:,1]))<0
 result={'guardMinJ':float(guardJ.min()),'guardMaxJ':float(guardJ.max()),'guardNonPositive':int((guardJ<=0).sum()),'sourceOpposedFaces':int(opposed.sum()),'experimental':True,'accepted':False,'target':target,'progress':progress,'bulk':bulk,'codeSha256':code_hashes,'inputSha256':hashlib.sha256((p/'cage.npz').read_bytes()).hexdigest(),'cageNodes':len(rest),'sourceVertices':len(source),'surfaceVolumeRatio':vol(now)/vol(source),'surfaceEdgeRatio':[float(er.min()),float(er.max())],'surfaceAreaRatio':[float((area1/area0).min()),float((area1/area0).max())],'attachmentMaxMm':error,'allStagesConverged':all(r['converged']for r in stages),'seconds':time.time()-clock,'stages':stages,'boneContactEvaluated':False,'selfCollisionEvaluated':False,'activation':None,'force':None}
 out.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');np.savez_compressed(out.with_suffix('.npz'),nodes=x,referenceNodes=rest,spacingMm=a['spacingMm'],origin=a['origin'],shape=a['shape'],sourceWorld=(now*scale+center)@RA.T+ta,referenceWorld=source*scale+center,sourceFaces=faces,sourceWeights=a['sourceWeights'])
 if not report['converged']:raise RuntimeError('Cage stage did not converge; partial result saved')
 step=min(base_step,step*1.25)
print('FINISHED',json.dumps({k:v for k,v in result.items()if k!='stages'}),flush=True)
