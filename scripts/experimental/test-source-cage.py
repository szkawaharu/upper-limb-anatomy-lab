import sys,json,tempfile,subprocess
from pathlib import Path
import numpy as np
with tempfile.TemporaryDirectory()as temporary:
 root=Path(temporary);src=root/'input';out=root/'output';src.mkdir();p=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1],[.25,.25,.25]]);nodes=np.r_[p,[[9,9,9]]];t=np.array([[4,1,2,3],[0,4,2,3],[0,1,4,3],[0,1,2,4]]);f=np.array([[1,2,3],[0,3,2],[0,1,3],[0,2,1]])
 np.savez(src/'fixture-source.npz',p=p,faces=f,weights=[0,0,1,1,.5]);np.savez(src/'fixture-tetrahedra.npz',nodes=nodes,tets=t);(src/'targets.json').write_text(json.dumps({'id':'fixture','worldUnitsPerMm':1,'poses':[]}))
 subprocess.run([sys.executable,str(Path(__file__).with_name('prepare-source-cage.py')),str(src),str(out)],check=True,capture_output=True)
 a=np.load(out/'cage.npz');s=np.load(out/'conforming.npz');assert np.array_equal(a['tets'][a['embeddingElements']],a['embeddingVertices'])
 assert np.max(np.abs(np.einsum('ni,nij->nj',a['embeddingBarycentric'],a['nodes'][a['embeddingVertices']])-a['sourcePoints']))<1e-15
 for face,element in zip(f,s['elements']):assert set(face).issubset(set(t[element]))
 assert np.array_equal(s['classes'][:5],[0,0,1,1,-1])and len(a['nodes'])==6
print('PASS full-source adapter: correct incident element per original vertex, exact identity embedding, correct boundary-face ownership, inactive generated node retained')
