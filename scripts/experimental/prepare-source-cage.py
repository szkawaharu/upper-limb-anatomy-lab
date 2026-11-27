"""Identity embedding of an audited full source solid for shared solver tests.
No coarse interpolation or new anatomical attachment region is introduced.
This is an offline reference and may be far too large for real-time use.
"""
import sys,json,hashlib,shutil
from pathlib import Path
import numpy as np
src,out=map(Path,sys.argv[1:3]);out.mkdir(parents=True,exist_ok=True);targets=json.loads((src/'targets.json').read_text());prefix='deltoid-mid'if targets['id']=='v60DeltMid'else targets['id'];sf=src/(prefix+'-source.npz');mf=src/(prefix+'-tetrahedra.npz');s=np.load(sf);m=np.load(mf);scale=targets['worldUnitsPerMm'];center=s['p'].mean(0);nodes=(m['nodes']-center)/scale;points=(s['p']-center)/scale;n=len(points);faces=s['faces'];assert np.array_equal(m['nodes'][:n],s['p']);classes=np.full(len(nodes),-1,np.int8);fixed=np.flatnonzero(np.isin(s['weights'],[0,1]));classes[fixed]=s['weights'][fixed];tets=m['tets'];Dm=np.stack([nodes[tets[:,k]]-nodes[tets[:,0]]for k in [1,2,3]],axis=2);volume=np.linalg.det(Dm)/6;assert np.all(volume>0);owner={};incident=np.full(len(nodes),-1,int)
for i,t in enumerate(tets):
 incident[t]=i
 for j in range(4):owner[tuple(sorted(np.delete(t,j)))]=i
assert np.all(incident[:n]>=0);embedding=tets[incident[:n]];bary=(embedding==np.arange(n)[:,None]).astype(float);assert np.all(bary.sum(1)==1)
elements=np.array([owner[tuple(sorted(f))]for f in faces]);np.savez_compressed(out/'cage.npz',nodes=nodes,tets=tets,weights=volume/volume.sum(),classes=classes,sourcePoints=points,sourceFaces=faces,sourceWeights=s['weights'],embeddingVertices=embedding,embeddingBarycentric=bary,embeddingElements=incident[:n],centerWorld=center,worldUnitsPerMm=scale)
np.savez_compressed(out/'conforming.npz',points=points,faces=faces,weights=s['weights'],classes=classes,embeddingVertices=embedding,embeddingBarycentric=bary,sourceTriangle=np.arange(len(faces)),elements=elements)
report={'experimental':True,'accepted':False,'id':targets['id'],'method':'full source tetrahedra with identity surface embedding; not real-time','sourceVertices':n,'sourceTriangles':len(faces),'nodes':len(nodes),'tets':len(tets),'fixedNodes':len(fixed),'sourceOriginalVerticesUnchanged':True,'sourceAttachmentMasksUnchanged':True,'volumeIntegration':'exact tetrahedral element volume','inputs':{str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in [sf,mf]}};(out/'preparation.json').write_text(json.dumps(report,indent=2)+'\n');shutil.copy2(src/'targets.json',out/'targets.json');print(json.dumps(report))
