"""Stateless positive-J continuation for the reduced source solid. This is a
numerical boundary-loading path, not a claim about measured joint motion.
"""
import sys,json,time,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from reduced_volume_model import ReducedSolid,affine_boundary_lift
directory=Path(sys.argv[1]);count=int(sys.argv[2]);elev=float(sys.argv[3]);plane=float(sys.argv[4]);bulk=float(sys.argv[5]if len(sys.argv)>5 else 30);iterations=int(sys.argv[6]if len(sys.argv)>6 else 40)
file=directory/('reduced-'+str(count)+'.npz');data=np.load(file);targets=json.loads((directory/'targets.json').read_text());target=next(x for x in targets['poses']if x['elev']==elev and x['plane']==plane);scale=float(data['worldUnitsPerMm']);center=data['centerWorld'];reference=data['reference'];clock=time.time();model=ReducedSolid(reference,data['tets'],data['fixed'],data['basis'])
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
 A,at=affine(target['a'],angles[0],progress);B,bt=affine(target['b'],angles[1],progress);return affine_boundary_lift(reference,data['lifting'],A,at,B,bt)
progress=0.;base_step=min(1.,np.deg2rad(3)/max(max(x[1]for x in angles),1e-12));step=base_step;q=np.zeros((data['basis'].shape[1],3));rows=[]
prefix=directory/('reduced-cont-'+str(count)+'-e'+str(int(elev))+'-p'+str(int(plane))+'-k'+str(int(bulk)))
while progress<1-1e-12:
 next_progress=min(1.,progress+step);base=base_at(next_progress)
 if not np.isfinite(model.evaluate(base,q,bulk)[0]):
  step*=.5
  if step<1e-7:raise RuntimeError('Reduced space cannot advance positive-J boundary loading')
  continue
 x,q,r=model.minimize(base,q,bulk=bulk,iterations=iterations)
 progress=next_progress;rows.append({'progress':progress,'step':step,'elapsedSeconds':time.time()-clock,**r});print(json.dumps({k:v for k,v in rows[-1].items()if k!='rows'}),flush=True)
 result={'experimental':True,'accepted':False,'target':target,'progress':progress,'allStagesConverged':all(x['converged']for x in rows),'attachmentMaxMm':float(np.linalg.norm(x[data['fixed']]-base[data['fixed']],axis=1).max()),'seconds':time.time()-clock,'basisSha256':hashlib.sha256(file.read_bytes()).hexdigest(),'codeSha256':{p:hashlib.sha256(Path(__file__).with_name(p).read_bytes()).hexdigest()for p in ['reduced_volume_model.py','volumetric_reference.py','run-reduced-continuation.py']},'stages':rows}
 prefix.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');np.savez(prefix.with_suffix('.npz'),world=x*scale+center,referenceWorld=reference*scale+center,tets=data['tets'],sourceFaces=data['faces'],fixed=data['fixed'],q=q);step=min(base_step,step*1.25)
print('FINISHED',json.dumps({'progress':progress,'seconds':time.time()-clock,'allStagesConverged':result['allStagesConverged']}),flush=True)
