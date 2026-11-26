"""Pose-independent Dirichlet harmonic bone blend fields. Original exact
attachment masks remain unchanged. Components without BOTH anchor values keep
their prior estimated field, rather than inventing new anatomical attachment.
"""
import json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import cg
import sys
src=Path(sys.argv[1]if len(sys.argv)>1 else'test-results/v076-source-clean/bindings.json');out=Path(sys.argv[2]if len(sys.argv)>2 else'test-results/v076-source-harmonic');out.mkdir(parents=True,exist_ok=True)
records=[];bindings=[]
for B in json.loads(src.read_text()):
 p=np.array(B['rest']).reshape(-1,3);f=np.array(B['faces']).reshape(-1,3);old=np.array(B['weights']);n=len(old);edges=np.unique(np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1),axis=0);a,b=edges.T;value=1/np.maximum(np.linalg.norm(p[a]-p[b],axis=1),1e-6);W=coo_matrix((np.r_[value,value],(np.r_[a,b],np.r_[b,a])),shape=(n,n)).tocsr();L=diags(np.asarray(W.sum(1)).ravel())-W;count,labels=connected_components(W);result=old.copy();logs=[]
 for k in range(count):
  ids=np.flatnonzero(labels==k);fixed=ids[(old[ids]==0)|(old[ids]==1)];free=ids[(old[ids]>0)&(old[ids]<1)];values=np.unique(old[fixed])
  if len(values)<2:logs.append({'vertices':len(ids),'free':len(free),'status':'retained-insufficient-attachment-data'});continue
  if not len(free):continue
  A=L[free][:,free];rhs=-(L[free][:,fixed]@old[fixed]);x,code=cg(A,rhs,x0=old[free],M=diags(1/A.diagonal()),rtol=1e-10,atol=1e-12,maxiter=20000)
  if code!=0:raise RuntimeError((B['spec']['id'],code))
  result[free]=np.clip(x,2**-149,1-2**-24);residual=float(np.linalg.norm(A@x-rhs)/max(np.linalg.norm(rhs),1e-20));logs.append({'vertices':len(ids),'free':len(free),'status':'harmonic','relativeResidual':residual})
 result=result.astype('<f4');assert np.array_equal(result==0,old==0)and np.array_equal(result==1,old==1),'Dirichlet masks changed'
 B['weights']=result.tolist();B['material']=result.tolist();B['harmonicAttachments']={'method':'positive inverse-edge graph Dirichlet harmonic solve; no fidelity term','components':logs,'fixedVertices':int(((old==0)|(old==1)).sum()),'maxRestFieldChange':float(np.max(np.abs(result-old))),'poseDependent':False};bindings.append(B);records.append({'id':B['spec']['id'],**B['harmonicAttachments']});print(json.dumps(records[-1]),flush=True)
(out/'bindings.json').write_text(json.dumps(bindings,separators=(',',':')));(out/'report.json').write_text(json.dumps({'method':__doc__,'meshes':records},indent=2)+'\n')
