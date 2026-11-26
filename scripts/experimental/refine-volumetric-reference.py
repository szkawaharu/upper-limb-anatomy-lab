"""Refine one saved feasible pose at an explicitly recorded volume penalty.
This changes no Dirichlet boundary or geometry. Not a real-time implementation.
"""
import sys,json,time,hashlib
from pathlib import Path
import numpy as np
from volumetric_reference import Solid
input=Path(sys.argv[1]);bulk=float(sys.argv[2]);iterations=int(sys.argv[3]if len(sys.argv)>3 else 500);out=Path(sys.argv[4])
a=np.load(input);scale=.00555;center=a['referenceWorld'].mean(0);rest=(a['referenceWorld']-center)/scale;start=(a['world']-center)/scale;solid=Solid(rest,a['tets'],a['fixed']);clock=time.time()
x,r=solid.minimize(start,bulk=bulk,iterations=iterations)
r.update({'input':str(input),'inputSha256':hashlib.sha256(input.read_bytes()).hexdigest(),'seconds':time.time()-clock,'fixedDeltaMm':float(np.linalg.norm(x[solid.fixed]-start[solid.fixed],axis=1).max()),'codeSha256':hashlib.sha256(Path(__file__).with_name('volumetric_reference.py').read_bytes()).hexdigest(),'experimental':True,'accepted':False})
np.savez(out.with_suffix('.npz'),world=x*scale+center,tets=a['tets'],referenceWorld=a['referenceWorld'],sourceFaces=a['sourceFaces'],fixed=a['fixed']);out.with_suffix('.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
