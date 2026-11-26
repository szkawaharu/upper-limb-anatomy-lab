"""Offline source-solid response basis under numerical boundary continuation.
Results always remain experimental; a reduced gradient cannot certify anatomy.
"""
import sys,json,time,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from reduced_volume_model import affine_boundary_lift
from vector_reduced_solid import VectorReducedSolid

directory=Path(sys.argv[1]);basis_name=sys.argv[2];elev=float(sys.argv[3]);plane=float(sys.argv[4]);bulk=float(sys.argv[5]if len(sys.argv)>5 else 30);iterations=int(sys.argv[6]if len(sys.argv)>6 else 40)
file=directory/(basis_name+'.npz');data=np.load(file);targets=json.loads((directory/'targets.json').read_text());target=next(x for x in targets['poses']if x['elev']==elev and x['plane']==plane);scale=float(data['worldUnitsPerMm']);center=data['centerWorld'];reference=data['reference'];clock=time.time();model=VectorReducedSolid(reference,data['tets'],data['fixed'],data['vectorBasis'])
code_hashes={p:hashlib.sha256(Path(__file__).with_name(p).read_bytes()).hexdigest()for p in ['vector_reduced_solid.py','volumetric_reference.py','reduced_volume_model.py','run-vector-continuation.py']};basis_hash=hashlib.sha256(file.read_bytes()).hexdigest()
def angle_axis(bone):
 q=np.array(bone['q']);s=np.linalg.norm(q[:3]);a=2*np.arctan2(s,q[3])
 if s<1e-12:
  if q[3]<0:raise ValueError('Ambiguous full-turn lift')
  return np.zeros(3),0.
 return q[:3]/s,float(a)
angles=[angle_axis(target[k])for k in['a','b']]
def affine(bone,aa,progress):
 axis,angle=aa;R=Rotation.from_rotvec(axis*angle*progress).as_matrix();rest=np.array(bone['rest']);at=rest+progress*(np.array(bone['p'])-rest);return R,(R@center+at-R@rest-center)/scale

def base_at(progress):
 A,at=affine(target['a'],angles[0],progress);B,bt=affine(target['b'],angles[1],progress)
 # Source-coordinate solve avoids artificial world rotation of modal vectors.
 local=affine_boundary_lift(reference,data['lifting'],np.eye(3),np.zeros(3),A.T@B,A.T@(bt-at))
 return local,A,at

progress=0.;base_step=min(1.,np.deg2rad(3)/max(max(x[1]for x in angles),1e-12));step=base_step;q=np.zeros(data['vectorBasis'].shape[2]);rows=[]
prefix=directory/('vector-cont-'+basis_name+'-e'+str(int(elev))+'-p'+str(int(plane))+'-k'+str(int(bulk)))
while progress<1-1e-12:
 next_progress=min(1.,progress+step);base,A,at=base_at(next_progress)
 if not np.isfinite(model.evaluate(base,q,bulk)[0]):
  step*=.5
  if step<1e-7:raise RuntimeError('Response space cannot advance positive-J boundary loading')
  continue
 x,q,r=model.minimize(base,q,bulk=bulk,iterations=iterations)
 progress=next_progress;rows.append({'progress':progress,'step':step,'elapsedSeconds':time.time()-clock,**r});print(json.dumps({k:v for k,v in rows[-1].items()if k!='rows'}),flush=True)
 world=(x@A.T+at)*scale+center
 result={'experimental':True,'accepted':False,'target':target,'progress':progress,'allReducedStagesConverged':all(x['converged']for x in rows),'allFullStagesConverged':all(x['fullReferenceConverged']for x in rows),'attachmentMaxMm':float(np.linalg.norm(x[data['fixed']]-base[data['fixed']],axis=1).max()),'seconds':time.time()-clock,'basisSha256':basis_hash,'codeSha256':code_hashes,'stages':rows}
 prefix.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');np.savez(prefix.with_suffix('.npz'),world=world,referenceWorld=reference*scale+center,tets=data['tets'],sourceFaces=data['faces'],fixed=data['fixed'],q=q);step=min(base_step,step*1.25)
print('FINISHED',json.dumps({'progress':progress,'seconds':time.time()-clock,'allReducedStagesConverged':result['allReducedStagesConverged'],'allFullStagesConverged':result['allFullStagesConverged']}),flush=True)
