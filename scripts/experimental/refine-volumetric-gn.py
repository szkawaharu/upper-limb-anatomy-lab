import sys,json,time,hashlib
from pathlib import Path
import numpy as np
from volumetric_reference import Solid
from volumetric_gauss_newton import minimize
from volumetric_preconditioner import factor_rest_metric,factor_strain_metric,warped_rest_preconditioner
from scipy.sparse.linalg import LinearOperator
input=Path(sys.argv[1]);bulk=float(sys.argv[2]);iterations=int(sys.argv[3]);out=Path(sys.argv[4]);a=np.load(input);scale=.00555;center=a['referenceWorld'].mean(0);rest=(a['referenceWorld']-center)/scale;start=(a['world']-center)/scale;solid=Solid(rest,a['tets'],a['fixed']);clock=time.time()
code_hashes={p:hashlib.sha256(Path(__file__).with_name(p).read_bytes()).hexdigest()for p in ['volumetric_reference.py','volumetric_gauss_newton.py','volumetric_preconditioner.py','refine-volumetric-gn.py']}
preconditioner=None;preconditioner_info=None;factory=None
if len(sys.argv)>5 and sys.argv[5]in ['coupled','warped','frozen']:
 factor,preconditioner_info=factor_strain_metric(solid,start,bulk)if sys.argv[5]=='frozen'else factor_rest_metric(solid,bulk);n=len(solid.free)*3;preconditioner=LinearOperator((n,n),matvec=factor.solve,dtype=float);preconditioner_info['setupSeconds']=time.time()-clock;preconditioner_info['pose']='initial-current'if sys.argv[5]=='frozen'else'source-rest';print('PRECONDITIONER',json.dumps(preconditioner_info),flush=True)
 if sys.argv[5]=='warped':
  factory=lambda x:warped_rest_preconditioner(solid,x,factor);preconditioner_info['transport']='source-volume-averaged proper nodal polar rotations; energy unchanged'
def progress(row):print(json.dumps({'seconds':time.time()-clock,**row}),flush=True)
x,r=minimize(solid,start,bulk=bulk,iterations=iterations,callback=progress,preconditioner=preconditioner,preconditioner_info=preconditioner_info,preconditioner_factory=factory)
r.update({'input':str(input),'inputSha256':hashlib.sha256(input.read_bytes()).hexdigest(),'seconds':time.time()-clock,'fixedDeltaMm':float(np.linalg.norm(x[solid.fixed]-start[solid.fixed],axis=1).max()),'codeSha256':code_hashes,'experimental':True,'accepted':False})
np.savez(out.with_suffix('.npz'),world=x*scale+center,tets=a['tets'],referenceWorld=a['referenceWorld'],sourceFaces=a['sourceFaces'],fixed=a['fixed']);out.with_suffix('.json').write_text(json.dumps(r,indent=2)+'\n');print('FINISHED',json.dumps(r),flush=True)
