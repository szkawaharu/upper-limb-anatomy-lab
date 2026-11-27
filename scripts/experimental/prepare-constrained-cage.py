# SPDX-License-Identifier: AGPL-3.0-or-later
"""OFFLINE TetGen PLC generator. Source attachment points AND whole fixed
facets are recovered inside a coarse outer box. Neither original tissue nor
its fixed footprint is decimated. AGPL TetGen is not linked to the web app.
"""
import sys,json,time,hashlib,shutil,itertools
from pathlib import Path
import numpy as np,tetgen
from constrained_cage import ConstrainedCage
from conforming_surface import conform
from attachment_refinement import refine
src=Path(sys.argv[1]);spacing=float(sys.argv[2]);out=Path(sys.argv[3]);out.mkdir(parents=True,exist_ok=True);clock=time.time();surface=np.load(src/'deltoid-mid-source.npz');material=np.load(src/'deltoid-mid-tetrahedra.npz');scale=.00555;center=surface['p'].mean(0);points=(surface['p']-center)/scale;w=surface['weights'];fixed=np.flatnonzero(np.isin(w,[0,1]));mapping=np.full(len(w),-1,int);mapping[fixed]=np.arange(len(fixed))
low=np.floor(points.min(0)/spacing)*spacing-spacing;high=np.ceil(points.max(0)/spacing)*spacing+spacing;grid=np.stack(np.meshgrid(*[np.arange(low[k]+spacing,high[k],spacing)for k in range(3)],indexing='ij'),axis=-1).reshape(-1,3);corners=np.array(list(itertools.product(*zip(low,high))));nodes=np.r_[points[fixed],grid,corners];start=len(nodes)-8
box=np.array([[0,1,3],[0,3,2],[4,6,7],[4,7,5],[0,4,5],[0,5,1],[2,3,7],[2,7,6],[0,2,6],[0,6,4],[1,5,7],[1,7,3]])+start
F=surface['faces'];fixed_faces=F[np.all(w[F]==0,axis=1)|np.all(w[F]==1,axis=1)];facets=np.r_[box,mapping[fixed_faces]]
generator=tetgen.TetGen(nodes,facets);result=generator.tetrahedralize(plc=True,convex=True,quality=True,minratio=2.,steinerleft=20000,nobisect=True,nomergevertex=True,nojettison=True,docheck=True,quiet=True,facesout=True)
generated,tets=result[:2];assert np.array_equal(generated[:len(nodes)],nodes),'Input control points changed';classes=np.full(len(generated),-1,np.int8);classes[:len(fixed)]=w[fixed]
all_faces=np.sort(np.concatenate([tets[:,[1,2,3]],tets[:,[0,3,2]],tets[:,[0,1,3]],tets[:,[0,2,1]]]),axis=1);keys=set(map(tuple,all_faces));assert all(tuple(sorted(t))in keys for t in mapping[fixed_faces]),'An attachment facet was not recovered exactly'
source_edges=np.unique(np.sort(np.concatenate([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]]),axis=1),axis=0);fixed_edges=source_edges[(w[source_edges[:,0]]==w[source_edges[:,1]])&np.isin(w[source_edges[:,0]],[0,1])]
if len(sys.argv)>4 and sys.argv[4]=='full-refinement':generated,tets,classes,refinement,_=refine(generated,tets,classes,mapping[fixed_faces],mapping[fixed_edges])
else:refinement={'status':'deferred to source-boundary-directed refinement','poseIndependent':True}
c=ConstrainedCage(generated,tets,classes);classes,vertices,bary,ids=c.fixed_patches(points,w);weights,info=c.material_weights((material['nodes']-center)/scale,material['tets'])
np.savez_compressed(out/'cage.npz',nodes=c.nodes,tets=c.tets,weights=weights,classes=classes,sourcePoints=points,sourceFaces=F,sourceWeights=w,embeddingVertices=vertices,embeddingBarycentric=bary,embeddingElements=ids,centerWorld=center,worldUnitsPerMm=scale)
report={'experimental':True,'accepted':False,'method':'PLC cage with exact original fixed points and attachment facets','spacingMm':spacing,'nodes':len(c.nodes),'tets':len(c.tets),'sourceVertices':len(points),'sourceTriangles':len(F),'fixedNodes':int((classes>=0).sum()),'sourceAttachmentVertices':len(fixed),'fixedSourceFacets':len(fixed_faces),'fixedFacetsExact':True,'tetgenWrapperVersion':tetgen.__version__,'generatorLicense':'AGPL-3.0-or-later','runtimeDependencyAdded':False,'materialQuadrature':info,'sourceDirectory':str(src),'inputs':{str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in [src/'deltoid-mid-source.npz',src/'deltoid-mid-tetrahedra.npz']}}
report['attachmentRefinement']=refinement;(out/'preparation.json').write_text(json.dumps(report,indent=2)+'\n');shutil.copy2(src/'targets.json',out/'targets.json');print(json.dumps(report),flush=True)
r=conform(c,points,F,w,lambda i:print('face',i,flush=True));classes,ids,bary,t=c.fixed_patches(r['points'],r['weights']);provenance=r['newVertexProvenance'];np.savez_compressed(out/'conforming.npz',points=r['points'],faces=r['faces'],weights=r['weights'],classes=classes,embeddingVertices=ids,embeddingBarycentric=bary,sourceTriangle=r['sourceTriangle'],elements=r['element'],newVertexSourceTriangle=np.array([x[1]for x in provenance]),newVertexSourceBarycentric=np.array([x[2]for x in provenance]))
edges=np.sort(np.concatenate([r['faces'][:,[0,1]],r['faces'][:,[1,2]],r['faces'][:,[2,0]]]),axis=1);_,count=np.unique(edges,axis=0,return_counts=True);assert np.all(count==2),'Subdivision became nonmanifold'
report={k:r[k]for k in ['sourceAreaMaxRelativeError','maxMergeDistanceMm']};report.update({'experimental':True,'accepted':False,'sourceVertices':len(points),'conformingVertices':len(r['points']),'sourceTriangles':len(F),'conformingTriangles':len(r['faces']),'edgeMultiplicities':{str(k):int((count==k).sum())for k in np.unique(count)},'seconds':time.time()-clock});(out/'conforming.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
