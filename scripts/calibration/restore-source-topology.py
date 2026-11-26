"""Restore source faces and remove exact zero-area/zero-thickness artifacts.
The old parser kept only one face from coincident opposite pairs, opening a
closed source chain. Cancel both opposite faces, retain every other triangle,
then remove unreferenced vertices. No hole filling or anatomical remodelling.
"""
import json,struct,sys,hashlib
from pathlib import Path
import numpy as np,trimesh
src=Path('test-results/v076-source-bindings/bindings.json');out=Path('test-results/v076-source-clean');out.mkdir(parents=True,exist_ok=True)
records=[];bindings=[];cache=Path('/tmp/upper-limb-atlas')
for B in json.loads(src.read_text()):
 raw=(cache/B['spec']['file']).read_bytes();n=struct.unpack_from('<I',raw,80)[0];dtype=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')]);points=np.frombuffer(raw,offset=84,count=n,dtype=dtype)['v'].reshape(-1,3)
 unique,first,inverse=np.unique(points,axis=0,return_index=True,return_inverse=True);order=np.argsort(first);new=np.empty(len(order),int);new[order]=np.arange(len(order));F=new[inverse].reshape(-1,3)
 assert len(order)==len(B['weights'])
 groups={};remove=set();zero=0
 for i,tri in enumerate(F):
  if len(set(tri))<3:remove.add(i);zero+=1;continue
  key=tuple(sorted(tri));groups.setdefault(key,[]).append(i)
 pairs=0
 for key,ids in groups.items():
  if len(ids)!=2:continue
  a,b=F[ids];parity=lambda t:sum(t[i]>t[j]for i in range(3)for j in range(i+1,3))%2
  if parity(a)!=parity(b):remove.update(ids);pairs+=1
 keep=np.array([i for i in range(len(F))if i not in remove]);F=F[keep]
 used=np.unique(F);remap=np.full(len(order),-1,int);remap[used]=np.arange(len(used));F=remap[F];p=np.array(B['rest']).reshape(-1,3)[used]
 before=trimesh.Trimesh(vertices=np.array(B['rest']).reshape(-1,3),faces=np.array(B['faces']).reshape(-1,3),process=False)
 m=trimesh.Trimesh(vertices=p,faces=F,process=False)
 B['rest']=p.ravel().tolist();B['faces']=F.ravel().tolist()
 for key in ['weights','lane','material']:B[key]=np.array(B[key])[used].tolist()
 B['sourceVertexIds']=used.tolist()
 B['meshPreparation']={'method':'exact source weld; restore original triangle chain; cancel coincident opposite faces; discard zero-area index faces and unreferenced vertices','sourceSha256':hashlib.sha256(raw).hexdigest(),'sourceTriangles':n,'cancelledOppositeFacePairs':pairs,'zeroIndexTriangles':zero,'removedUnreferencedVertices':len(order)-len(used),'retainedTriangles':len(F),'closed':bool(m.is_watertight),'oriented':bool(m.is_winding_consistent)}
 records.append({'id':B['spec']['id'],'beforeClosed':bool(before.is_watertight),**B['meshPreparation']});bindings.append(B);print(json.dumps(records[-1]),flush=True)
(out/'bindings.json').write_text(json.dumps(bindings,separators=(',',':')))
(out/'report.json').write_text(json.dumps({'method':__doc__,'inputSha256':hashlib.sha256(src.read_bytes()).hexdigest(),'meshes':records},indent=2)+'\n')
