import sys,json,time,hashlib,shutil
from pathlib import Path
import numpy as np
from point_cage import PointCage
from conforming_surface import conform

src=Path(sys.argv[1]);spacing=float(sys.argv[2]);out=Path(sys.argv[3]);out.mkdir(parents=True,exist_ok=True);clock=time.time();surface=np.load(src/'deltoid-mid-source.npz');material=np.load(src/'deltoid-mid-tetrahedra.npz');scale=.00555;center=surface['p'].mean(0);points=(surface['p']-center)/scale
c=PointCage(points,surface['weights'],spacing);classes,vertices,bary,ids=c.fixed_patches(points,surface['weights']);weights,info=c.material_weights((material['nodes']-center)/scale,material['tets'])
np.savez_compressed(out/'cage.npz',nodes=c.nodes,tets=c.tets,weights=weights,classes=classes,sourcePoints=points,sourceFaces=surface['faces'],sourceWeights=surface['weights'],embeddingVertices=vertices,embeddingBarycentric=bary,embeddingElements=ids,centerWorld=center,worldUnitsPerMm=scale)
report={'experimental':True,'accepted':False,'method':'exact source attachment controls; unstructured Delaunay background','spacingMm':spacing,'nodes':len(c.nodes),'tets':len(c.tets),'sourceVertices':len(points),'sourceTriangles':len(surface['faces']),'fixedNodes':int((classes>=0).sum()),'sourceAttachmentVertices':int(np.isin(surface['weights'],[0,1]).sum()),'materialQuadrature':info,'sourceDirectory':str(src),'inputs':{str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in [src/'deltoid-mid-source.npz',src/'deltoid-mid-tetrahedra.npz']}}
(out/'preparation.json').write_text(json.dumps(report,indent=2)+'\n');shutil.copy2(src/'targets.json',out/'targets.json');print(json.dumps(report),flush=True)
r=conform(c,points,surface['faces'],surface['weights'],lambda i:print('face',i,flush=True))
try:classes,ids,bary,t=c.fixed_patches(r['points'],r['weights'])
except ValueError as e:
 (out/'rejected.json').write_text(json.dumps({'experimental':True,'accepted':False,'reason':str(e),'conclusion':'Unconstrained point Delaunay does not preserve the whole declared attachment area','sourceAreaMaxRelativeError':r['sourceAreaMaxRelativeError']},indent=2)+'\n');raise
provenance=r['newVertexProvenance'];np.savez_compressed(out/'conforming.npz',points=r['points'],faces=r['faces'],weights=r['weights'],classes=classes,embeddingVertices=ids,embeddingBarycentric=bary,sourceTriangle=r['sourceTriangle'],elements=r['element'],newVertexSourceTriangle=np.array([x[1]for x in provenance]),newVertexSourceBarycentric=np.array([x[2]for x in provenance]))
edges=np.sort(np.concatenate([r['faces'][:,[0,1]],r['faces'][:,[1,2]],r['faces'][:,[2,0]]]),axis=1);_,count=np.unique(edges,axis=0,return_counts=True)
report={k:r[k]for k in ['sourceAreaMaxRelativeError','maxMergeDistanceMm']};report.update({'experimental':True,'accepted':False,'sourceVertices':len(points),'conformingVertices':len(r['points']),'sourceTriangles':len(surface['faces']),'conformingTriangles':len(r['faces']),'edgeMultiplicities':{str(k):int((count==k).sum())for k in np.unique(count)},'seconds':time.time()-clock});(out/'conforming.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
