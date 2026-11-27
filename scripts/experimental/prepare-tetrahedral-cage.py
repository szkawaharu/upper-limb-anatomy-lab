import sys,json,hashlib
from pathlib import Path
import numpy as np
from tetrahedral_cage import GridCage
src=Path(sys.argv[1]);spacing=float(sys.argv[2]);out=Path(sys.argv[3]);out.mkdir(parents=True,exist_ok=True)
surface=np.load(src/'deltoid-mid-source.npz');material=np.load(src/'deltoid-mid-tetrahedra.npz');scale=.00555;center=surface['p'].mean(0);points=(surface['p']-center)/scale
cage=GridCage(points,spacing);classes,vertices,bary,ids=cage.fixed_patches(points,surface['weights']);weights,info=cage.material_weights((material['nodes']-center)/scale,material['tets'])
np.savez_compressed(out/'cage.npz',nodes=cage.nodes,tets=cage.tets,weights=weights,classes=classes,sourcePoints=points,sourceFaces=surface['faces'],sourceWeights=surface['weights'],embeddingVertices=vertices,embeddingBarycentric=bary,embeddingElements=ids,centerWorld=center,worldUnitsPerMm=scale)
report={'experimental':True,'accepted':False,'spacingMm':spacing,'nodes':len(cage.nodes),'tets':len(cage.tets),'sourceVertices':len(points),'sourceTriangles':len(surface['faces']),'fixedNodes':int((classes>=0).sum()),'sourceAttachmentVertices':int(np.isin(surface['weights'],[0,1]).sum()),'sourceReconstructionMaxMm':float(np.max(np.linalg.norm(cage.positions(cage.nodes,vertices,bary)-points,axis=1))),'materialQuadrature':info,'sourceDirectory':str(src),'inputs':{str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in [src/'deltoid-mid-source.npz',src/'deltoid-mid-tetrahedra.npz']}}
(out/'preparation.json').write_text(json.dumps(report,indent=2)+'\n');(out/'targets.json').write_bytes((src/'targets.json').read_bytes());print(json.dumps(report))
