"""Pose-independent, source-boundary-directed cage refinement.
Adds free controls only where undeclared source tissue has entirely fixed
interpolation supports. Inherits the exact parent-cell integration density.
No contact result, target angle, muscle-specific factor or score is consulted.
"""
import sys,json,time,hashlib,shutil
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from attachment_refinement import refine
from constrained_cage import ConstrainedCage
from conforming_surface import conform

src=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True);clock=time.time();a=np.load(src/'cage.npz');s=np.load(src/'conforming.npz');classes=s['classes'];ids=s['embeddingVertices'];bary=s['embeddingBarycentric'];weights=s['weights'];active=bary>1e-9;support_classes=classes[ids];unlabelled=~np.isin(weights,[0,1]);rigid=unlabelled&np.all((support_classes>=0)|~active,axis=1);needed={'edges':set(),'faces':set(),'cells':set()}
for i in np.flatnonzero(rigid):
 support=ids[i][active[i]];key=tuple(sorted(support))
 if len(support)==2:needed['edges'].add(key)
 elif len(support)==3:needed['faces'].add(key)
 elif len(support)==4:
  cell=np.flatnonzero(np.all(np.sort(a['tets'],axis=1)==np.sort(support),axis=1));needed['cells'].update(cell.tolist())
cl=classes[a['tets']];all_fixed=np.all(cl>=0,axis=1);source_cells=s['elements'][unlabelled[s['faces']].any(1)];needed['cells'].update(source_cells[all_fixed[source_cells]].tolist())
P=a['sourcePoints'];F=a['sourceFaces'];W=a['sourceWeights'];fixed_faces=F[np.all(W[F]==0,axis=1)|np.all(W[F]==1,axis=1)];e=np.unique(np.sort(np.concatenate([F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]]),axis=1),axis=0);fixed_edges=e[(W[e[:,0]]==W[e[:,1]])&np.isin(W[e[:,0]],[0,1])];distance,mapping=cKDTree(a['nodes']).query(P);assert distance[np.isin(W,[0,1])].max()==0
nodes,tets,classes,details,parents=refine(a['nodes'],a['tets'],classes,mapping[fixed_faces],mapping[fixed_edges],needed)
V=lambda p,t:np.linalg.det(np.stack([p[t[:,k]]-p[t[:,0]]for k in [1,2,3]],axis=2))/6
volume=V(nodes,tets);parent_volume=V(a['nodes'],a['tets']);integration=a['weights'][parents]*volume/parent_volume[parents];assert abs(integration.sum()-a['weights'].sum())<1e-12
c=ConstrainedCage(nodes,tets,classes);cl,emb,b,el=c.fixed_patches(P,W);np.savez_compressed(out/'cage.npz',**{**{k:a[k]for k in a.files},'nodes':nodes,'tets':tets,'weights':integration,'classes':cl,'embeddingVertices':emb,'embeddingBarycentric':b,'embeddingElements':el})
report=json.loads((src/'preparation.json').read_text());report.update({'method':'PLC cage with source-boundary-directed free-control refinement','nodes':len(nodes),'tets':len(tets),'fixedNodes':int((cl>=0).sum()),'attachmentRefinement':details,'unlabelledRigidSourceSamplesBefore':int(rigid.sum()),'targetedEdges':len(needed['edges']),'targetedFaces':len(needed['faces']),'targetedCells':len(needed['cells']),'parentCellDensityPreserved':True,'parentInputs':{name:hashlib.sha256((src/name).read_bytes()).hexdigest()for name in ['cage.npz','conforming.npz','preparation.json']}});(out/'preparation.json').write_text(json.dumps(report,indent=2)+'\n');shutil.copy2(src/'targets.json',out/'targets.json');print(json.dumps(report),flush=True)
r=conform(c,P,F,W,lambda i:print('face',i,flush=True));cl,emb,b,el=c.fixed_patches(r['points'],r['weights']);provenance=r['newVertexProvenance'];np.savez_compressed(out/'conforming.npz',points=r['points'],faces=r['faces'],weights=r['weights'],classes=cl,embeddingVertices=emb,embeddingBarycentric=b,sourceTriangle=r['sourceTriangle'],elements=r['element'],newVertexSourceTriangle=np.array([x[1]for x in provenance]),newVertexSourceBarycentric=np.array([x[2]for x in provenance]))
edges=np.sort(np.concatenate([r['faces'][:,[0,1]],r['faces'][:,[1,2]],r['faces'][:,[2,0]]]),axis=1);_,count=np.unique(edges,axis=0,return_counts=True);assert np.all(count==2)
report={k:r[k]for k in ['sourceAreaMaxRelativeError','maxMergeDistanceMm']};report.update({'experimental':True,'accepted':False,'sourceVertices':len(P),'conformingVertices':len(r['points']),'sourceTriangles':len(F),'conformingTriangles':len(r['faces']),'edgeMultiplicities':{str(k):int((count==k).sum())for k in np.unique(count)},'seconds':time.time()-clock});(out/'conforming.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
