import sys,json,time,hashlib
from pathlib import Path
import numpy as np
from volumetric_reference import Solid
from linear_elastic_modes import build_modes
directory=Path(sys.argv[1]);bulk=float(sys.argv[2]if len(sys.argv)>2 else 30);source=np.load(directory/'deltoid-mid-source.npz');data=np.load(directory/'reduced-24.npz');solid=Solid(data['reference'],data['tets'],data['fixed']);classes=source['weights'][solid.fixed].astype(int);clock=time.time()
def progress(row):print(json.dumps({'seconds':time.time()-clock,**row}),flush=True)
r=build_modes(solid,data['lifting'],classes,bulk,progress);prefix=directory/('elastic-modes-k'+str(int(bulk)));np.savez_compressed(prefix.with_suffix('.npz'),**{k:v for k,v in r.items()if k!='reports'},reference=data['reference'],tets=data['tets'],fixed=data['fixed'],faces=data['faces'],lifting=data['lifting'],centerWorld=data['centerWorld'],worldUnitsPerMm=data['worldUnitsPerMm'])
report={'experimental':True,'accepted':False,'bulk':bulk,'modes':r['vectorBasis'].shape[2],'discardedNumericalModes':r['discardedNumericalModes'],'reports':r['reports'],'seconds':time.time()-clock,'codeSha256':hashlib.sha256(Path(__file__).with_name('linear_elastic_modes.py').read_bytes()).hexdigest()};prefix.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps(report),flush=True)
