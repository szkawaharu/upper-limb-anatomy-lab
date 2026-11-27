"""Near-incompressible penalty continuation of a completed reference pose.
Only the numerical loading path changes. Each bulk value must reach the same
0.001mm equilibrium tolerance, with exact pins and positive-J Armijo checks.
Final constitutive parameters are not inferred muscle force/material data.
"""
import sys,json,time,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from volumetric_reference import Solid
from cage_solver import minimize
from snapshot_names import snapshot_name
file=Path(sys.argv[1]);r=json.loads(file.with_suffix('.json').read_text());a=np.load(file);assert r['progress']==1 and r['allStagesConverged'];assert hashlib.sha256(file.read_bytes()).hexdigest()==r['snapshotSha256'];p=file.parent;c=np.load(p/'cage.npz');s=np.load(p/'conforming.npz');fixed=np.flatnonzero(s['classes']>=0);solid=Solid(c['nodes'],c['tets'],fixed,element_weights=c['weights']);x=a['nodes'].copy();pins=x[fixed].copy();faces=s['faces'];source=s['points'];scale=float(c['worldUnitsPerMm']);center=c['centerWorld'];surface=lambda x:np.einsum('ni,nij->nj',s['embeddingBarycentric'],x[s['embeddingVertices']]);clock=time.time();history=[];RA=Rotation.from_quat(r['target']['a']['q']).as_matrix();ta=np.array(r['target']['a']['p'])-RA@np.array(r['target']['a']['rest']);bulk_values=[100.,300.,1000.,3000.];code={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()for name in ['refine-material-state.py','cage_solver.py','volumetric_reference.py','projected_solid_metric.py','constitutive_metric.py']}
for bulk in bulk_values:
 if bulk<=r['bulk']:continue
 output=p/snapshot_name(r['target']['elev'],r['target']['plane'],bulk)
 def observe(nodes,row):
  if row['iteration']%10==0:
   print(json.dumps({'bulk':bulk,'seconds':time.time()-clock,**row}),flush=True);temporary=p/'material-trial.tmp.npz';np.savez_compressed(temporary,nodes=nodes,bulk=bulk,converged=False);temporary.replace(p/'material-trial.npz')
 x,report=minimize(solid,x,bulk=bulk,iterations=60,observer=observe);assert np.array_equal(x[fixed],pins);history.append({'bulk':bulk,'seconds':time.time()-clock,**report});now=surface(x);R=source[faces];N=now[faces];n0=np.cross(R[:,1]-R[:,0],R[:,2]-R[:,0]);n=np.cross(N[:,1]-N[:,0],N[:,2]-N[:,0]);F=(solid.edges(x)@solid.inverse)[s['elements']];expected=np.linalg.det(F)[:,None]*np.einsum('tij,tj->ti',np.linalg.inv(F).transpose(0,2,1),n0);edges=np.unique(np.sort(np.r_[faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]],axis=1),axis=0);edge=np.linalg.norm(now[edges[:,1]]-now[edges[:,0]],axis=1)/np.linalg.norm(source[edges[:,1]]-source[edges[:,0]],axis=1);area=np.linalg.norm(n,axis=1)/np.linalg.norm(n0,axis=1);vol=lambda y:np.einsum('ti,ti->t',y[faces[:,0]],np.cross(y[faces[:,1]],y[faces[:,2]])).sum()/6
 result={**r,'bulk':bulk,'allStagesConverged':report['converged']and all(h['converged']for h in history),'materialContinuation':history,'materialTargetReached':bulk==bulk_values[-1]and report['converged'],'materialInitialization':{'input':str(file),'inputSha256':hashlib.sha256(file.read_bytes()).hexdigest(),'bulk':r['bulk']},'materialCodeSha256':code,'surfaceVolumeRatio':float(vol(now)/vol(source)),'surfaceEdgeRatio':[float(edge.min()),float(edge.max())],'surfaceAreaRatio':[float(area.min()),float(area.max())],'sourceOpposedFaces':int((np.einsum('ti,ti->t',n,expected)<0).sum()),'attachmentMaxMm':0.,'seconds':time.time()-clock,'stages':r['stages']+[{'progress':1.,'phase':'material-continuation',**h}for h in history]}
 temporary=output.with_suffix('.tmp.npz');np.savez_compressed(temporary,**{**{k:a[k]for k in a.files},'nodes':x,'sourceWorld':(now*scale+center)@RA.T+ta});result['snapshotSha256']=hashlib.sha256(temporary.read_bytes()).hexdigest();temporary.replace(output.with_suffix('.npz'));temporary=output.with_suffix('.tmp.json');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(output.with_suffix('.json'));print(json.dumps({'bulk':bulk,'converged':report['converged'],'residualMaxMm':report['residualMaxMm'],'JRange':[report['cageJQuantiles'][0],report['cageJQuantiles'][-1]],'seconds':time.time()-clock}),flush=True)
 if not report['converged']:raise RuntimeError('Material continuation did not converge; unaccepted checkpoint saved')
