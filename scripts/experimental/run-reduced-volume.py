import sys,json,time,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from reduced_volume_model import ReducedSolid
directory=Path(sys.argv[1]);count=int(sys.argv[2]);elev=float(sys.argv[3]);plane=float(sys.argv[4]);bulk=float(sys.argv[5]if len(sys.argv)>5 else 30);basisfile=directory/('reduced-'+str(count)+'.npz');data=np.load(basisfile);targets=json.loads((directory/'targets.json').read_text());target=next(x for x in targets['poses']if x['elev']==elev and x['plane']==plane);scale=float(data['worldUnitsPerMm']);center=data['centerWorld'];reference=data['reference'];world=reference*scale+center
def affine(bone):
 R=Rotation.from_quat(bone['q']).as_matrix();t=(R@center+np.array(bone['p'])-R@np.array(bone['rest'])-center)/scale;return R,t
A,at=affine(target['a']);B,bt=affine(target['b']);lifting=data['lifting']
base=reference@A.T+at+lifting[:,:3]@(B-A).T+lifting[:,3,None]*(bt-at)
clock=time.time();model=ReducedSolid(reference,data['tets'],data['fixed'],data['basis'])
def progress(row):print(json.dumps({'seconds':time.time()-clock,**row}),flush=True)
x,q,r=model.minimize(base,bulk=bulk,iterations=40,callback=progress)
r.update({'experimental':True,'accepted':False,'target':target,'seconds':time.time()-clock,'attachmentMaxMm':float(np.linalg.norm(x[data['fixed']]-base[data['fixed']],axis=1).max()),'basisSha256':hashlib.sha256(basisfile.read_bytes()).hexdigest(),'codeSha256':{p:hashlib.sha256(Path(__file__).with_name(p).read_bytes()).hexdigest()for p in ['reduced_volume_model.py','run-reduced-volume.py','volumetric_reference.py']}})
prefix=directory/('reduced-'+str(count)+'-e'+str(int(elev))+'-p'+str(int(plane))+'-k'+str(int(bulk)));np.savez(prefix.with_suffix('.npz'),world=x*scale+center,referenceWorld=world,tets=data['tets'],sourceFaces=data['faces'],fixed=data['fixed'],q=q);prefix.with_suffix('.json').write_text(json.dumps(r,indent=2)+'\n');print('FINISHED',json.dumps(r),flush=True)
