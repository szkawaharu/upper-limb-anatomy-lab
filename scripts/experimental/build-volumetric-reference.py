# SPDX-License-Identifier: AGPL-3.0-or-later
"""Optional OFFLINE generation tool, importing TetGen's AGPL core. This file
is separate from the web runtime and does not change other project licences.
Source anatomy retains BodyParts3D CC BY-SA 2.1 Japan attribution. No boundary
repair/decimation is attempted when source tetrahedralization fails.
"""
import hashlib,json,sys,time
from pathlib import Path
import numpy as np
import tetgen

input=Path(sys.argv[1]if len(sys.argv)>1 else'test-results/v076-source-refit/bindings.json')
id=sys.argv[2]if len(sys.argv)>2 else'v60DeltMid'
out=Path(sys.argv[3]if len(sys.argv)>3 else'test-results/v076-volumetric');out.mkdir(parents=True,exist_ok=True)
B=next(x for x in json.loads(input.read_text())if x['spec']['id']==id)
if not B['meshPreparation']['closed']or not B['meshPreparation']['oriented']:
 raise ValueError('Source must be closed and consistently oriented')
p=np.array(B['rest']).reshape(-1,3);faces=np.array(B['faces'],np.int32).reshape(-1,3)
prefix='deltoid-mid'if id=='v60DeltMid'else id
np.savez(out/(prefix+'-source.npz'),p=p,faces=faces,weights=np.array(B['weights']),sourceVertexIds=B['sourceVertexIds'])
started=time.time();generator=tetgen.TetGen(p,faces)
result=generator.tetrahedralize(plc=True,quality=False,nobisect=True,nomergevertex=True,nojettison=True,docheck=True,quiet=True,facesout=True)
nodes,cells=result[:2]
if not np.array_equal(nodes[:len(p)],p):raise ValueError('Original source nodes changed')
D=np.stack([nodes[cells[:,j]]-nodes[cells[:,0]]for j in[1,2,3]],axis=2);det=np.linalg.det(D)
if not np.all(det>0):raise ValueError('Non-positive reference tetrahedra')
all_faces=np.sort(np.concatenate([cells[:,[1,2,3]],cells[:,[0,3,2]],cells[:,[0,1,3]],cells[:,[0,2,1]]]),axis=1)
unique,count=np.unique(all_faces,axis=0,return_counts=True);boundary=unique[count==1]
order=lambda f:f[np.lexsort(f.T[::-1])]
if not np.array_equal(order(boundary),order(np.sort(faces,axis=1))):raise ValueError('Source boundary faces changed')
S=p[faces]-p.mean(0);surface=float(np.einsum('ij,ij->i',S[:,0],np.cross(S[:,1],S[:,2])).sum()/6)
volume=float(det.sum()/6);relative=abs(volume-surface)/abs(surface)
if relative>1e-10:raise ValueError('Source volume does not match tetrahedra')
np.savez(out/(prefix+'-tetrahedra.npz'),nodes=nodes,tets=cells)
report={'id':id,'inputSha256':hashlib.sha256(input.read_bytes()).hexdigest(),'source':B['definition']['source'],'sourceSha256':B['meshPreparation']['sourceSha256'],'tetgenWrapperVersion':tetgen.__version__,'runtimeDependencyAdded':False,'sourceVertices':len(p),'nodes':len(nodes),'tets':len(cells),'boundaryFaces':len(boundary),'sourceBoundaryExact':True,'sourcePrefixMaxError':0,'referenceTetSignCounts':{'positive':len(cells),'negative':0,'zero':0},'tetVolumeSum':volume,'sourceSurfaceVolume':surface,'relativeVolumeDifference':relative,'conditionNumberQuantiles':np.quantile(np.linalg.cond(D),[0,.5,.95,.99,1]).tolist(),'unreferencedGeneratedNodes':np.setdiff1d(np.arange(len(nodes)),np.unique(cells)).tolist(),'seconds':time.time()-started}
(out/'reference-mesh-audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
