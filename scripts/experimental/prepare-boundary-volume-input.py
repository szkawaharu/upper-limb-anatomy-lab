"""Reuse an audited reference solid with revised, versioned boundary DATA.
Every source point/face and every interior material tetrahedron is unchanged.
"""
import sys,json,hashlib,shutil
from pathlib import Path
import numpy as np
src,binding,out=map(Path,sys.argv[1:4]);out.mkdir(parents=True,exist_ok=True);a=np.load(src/'deltoid-mid-source.npz');B=next(b for b in json.loads(binding.read_text())if b['spec']['id']=='v60DeltMid')
assert np.array_equal(np.array(B['rest']).reshape(-1,3),a['p'])and np.array_equal(np.array(B['faces']).reshape(-1,3),a['faces'])
np.savez_compressed(out/'deltoid-mid-source.npz',**{k:np.array(B['weights'])if k=='weights'else a[k]for k in a.files})
for name in ['deltoid-mid-tetrahedra.npz','targets.json']:shutil.copy2(src/name,out/name)
(out/'boundary-provenance.json').write_text(json.dumps({'experimental':True,'accepted':False,'sourceGeometryChanged':False,'materialTetrahedraChanged':False,'inputs':{str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in [binding,src/'deltoid-mid-source.npz',src/'deltoid-mid-tetrahedra.npz']}},indent=2)+'\n')
