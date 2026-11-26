import sys,json,time,hashlib
from pathlib import Path
import numpy as np
from reduced_volume_basis import build_basis
directory=Path(sys.argv[1]if len(sys.argv)>1 else'test-results/v076-volumetric-quality');count=int(sys.argv[2]if len(sys.argv)>2 else 24);source=np.load(directory/'deltoid-mid-source.npz');mesh=np.load(directory/'deltoid-mid-tetrahedra.npz');scale=.00555;center=mesh['nodes'].mean(0);nodes=(mesh['nodes']-center)/scale;start=time.time();result=build_basis(nodes,mesh['tets'],source['weights'],count)
out=directory/('reduced-'+str(count));np.savez_compressed(out.with_suffix('.npz'),**result,reference=nodes,tets=mesh['tets'],faces=source['faces'],centerWorld=center,worldUnitsPerMm=scale)
report={'experimental':True,'productionEnabled':False,'sourceVertices':len(source['weights']),'nodes':len(nodes),'tets':len(mesh['tets']),'controls':len(result['controls']),'modes':result['basis'].shape[1],'patchVertices':result['controlPatchVertices'].tolist(),'discardedNumericalModes':result['discardedNumericalModes'],'sourceBoundaryExact':True,'fixedBasisMax':float(np.abs(result['basis'][result['fixed']]).max()),'seconds':time.time()-start,'codeSha256':hashlib.sha256(Path(__file__).with_name('reduced_volume_basis.py').read_bytes()).hexdigest()};out.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
