import sys,json,time
from pathlib import Path
import numpy as np
from tetrahedral_cage import GridCage
from conforming_surface import conform
p=Path(sys.argv[1]);a=np.load(p/'cage.npz');h=json.loads((p/'preparation.json').read_text())['spacingMm'];c=GridCage(a['sourcePoints'],h);assert np.array_equal(c.nodes,a['nodes']);clock=time.time();r=conform(c,a['sourcePoints'],a['sourceFaces'],a['sourceWeights'],lambda i:print('face',i,flush=True));classes,ids,bary,t=c.fixed_patches(r['points'],r['weights'])
np.savez_compressed(p/'conforming.npz',points=r['points'],faces=r['faces'],weights=r['weights'],classes=classes,embeddingVertices=ids,embeddingBarycentric=bary,sourceTriangle=r['sourceTriangle'],elements=r['element'])
edges=np.sort(np.concatenate([r['faces'][:,[0,1]],r['faces'][:,[1,2]],r['faces'][:,[2,0]]]),axis=1);_,count=np.unique(edges,axis=0,return_counts=True)
report={k:r[k]for k in ['sourceAreaMaxRelativeError','maxMergeDistanceMm']};report.update({'experimental':True,'accepted':False,'sourceVertices':len(a['sourcePoints']),'conformingVertices':len(r['points']),'sourceTriangles':len(a['sourceFaces']),'conformingTriangles':len(r['faces']),'edgeMultiplicities':{str(k):int((count==k).sum())for k in np.unique(count)},'seconds':time.time()-clock});(p/'conforming.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
