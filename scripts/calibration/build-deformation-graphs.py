"""Embedded deformation graphs from unchanged source surfaces. Node classes
separate fixed attachment regions from free tissue; no anatomical mesh is
replaced by the graph. Graph topology follows the source mesh components.
"""
import json,sys
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix,diags,csr_matrix
from scipy.sparse.csgraph import connected_components,reverse_cuthill_mckee
from scipy.sparse.linalg import cg
from scipy.spatial import cKDTree
src=Path(sys.argv[1]if len(sys.argv)>1 else'test-results/v076-source-refit/bindings.json');out=Path(sys.argv[2]if len(sys.argv)>2 else'test-results/v076-source-graphs');out.mkdir(parents=True,exist_ok=True)
records=[];bindings=[]
for B in json.loads(src.read_text()):
 id=B['spec']['id'];p=np.array(B['rest']).reshape(-1,3);f=np.array(B['faces']).reshape(-1,3);w=np.array(B['weights']);edges=np.unique(np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1),axis=0);length=np.linalg.norm(p[edges[:,0]]-p[edges[:,1]],axis=1);n=len(p)
 adjacency=coo_matrix((np.ones(len(edges)*2),(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(n,n)).tocsr();components,labels=connected_components(adjacency)
 if id=='v60PecAbd'and not(w==0).any():
  # This existing atlas part originates in abdominal aponeurotic tissue, not
  # a registered bone surface. Mark the source inferior terminal cap as a
  # LOW-CONFIDENCE trunk-frame boundary; do not call it measured enthesis.
  spacing=float(np.median(length));cap=p[:,1]<=p[:,1].min()+spacing;w[cap]=0;B['aponeuroticOriginEstimate']={'axis':'source inferior world Y','widthWorld':spacing,'vertices':int(cap.sum()),'confidence':'low; terminal cap only, aponeurotic footprint not acquired'}
  value=1/np.maximum(length,1e-6);W=coo_matrix((np.r_[value,value],(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(n,n)).tocsr();L=diags(np.asarray(W.sum(1)).ravel())-W
  fixed=np.flatnonzero((w==0)|(w==1));free=np.flatnonzero((w>0)&(w<1));A=L[free][:,free];rhs=-(L[free][:,fixed]@w[fixed]);x,code=cg(A,rhs,x0=w[free],M=diags(1/A.diagonal()),rtol=1e-10,atol=1e-12,maxiter=20000);assert code==0;w[free]=np.clip(x,2**-149,1-2**-24);w=w.astype('<f4').astype(float);B['weights']=w.tolist()
 fn=np.cross(p[f[:,1]]-p[f[:,0]],p[f[:,2]]-p[f[:,0]]);faceArea=np.linalg.norm(fn,axis=1)*.5;mass=np.zeros(n)
 for k in range(3):np.add.at(mass,f[:,k],faceArea/3)
 mass=np.maximum(mass,np.finfo(float).tiny);h=np.sqrt(faceArea.sum()/192);klass=np.where(w==0,0,np.where(w==1,1,2));keys=np.c_[labels,klass,np.floor((p-p.min(0))/h).astype(int)];unique,inverse=np.unique(keys,axis=0,return_inverse=True);m=len(unique);den=np.bincount(inverse,weights=mass,minlength=m);nodes=np.column_stack([np.bincount(inverse,weights=mass*p[:,k],minlength=m)/den for k in range(3)]);nodeW=np.bincount(inverse,weights=mass*w,minlength=m)/den;nodeClass=unique[:,1].copy();nodeComponent=unique[:,0]
 unanchored=[]
 for comp in range(components):
  ids=np.flatnonzero(nodeComponent==comp)
  if not np.isin(nodeClass[ids],[0,1]).any():
   nodeClass[ids]=3;nodeW[ids]=np.average(nodeW[ids],weights=den[ids]);unanchored.append({'component':int(comp),'nodes':len(ids),'vertices':int((labels==comp).sum()),'method':'rigid continuous bone blend; attachment data absent'})
 E=np.unique(np.sort(inverse[edges],axis=1),axis=0);E=E[E[:,0]!=E[:,1]];distance=np.linalg.norm(nodes[E[:,0]]-nodes[E[:,1]],axis=1);ew=1/np.maximum(distance,1e-6);W=coo_matrix((np.r_[ew,ew],(np.r_[E[:,0],E[:,1]],np.r_[E[:,1],E[:,0]])),shape=(m,m)).tocsr();L=diags(np.asarray(W.sum(1)).ravel())-W;free=np.flatnonzero(nodeClass==2);fixed=np.flatnonzero(nodeClass!=2)
 if len(free):
  A=L[free][:,free].tocsr();order=reverse_cuthill_mckee(A,symmetric_mode=True);free=free[order];A=L[free][:,free].toarray();chol=np.linalg.cholesky(A);lower=csr_matrix(np.tril(chol,-1));upper=lower.T.tocsr();diag=chol.diagonal()
 else:lower=upper=csr_matrix((0,0));diag=np.array([])
 embedIds=np.zeros((n,4),np.int32);embedWeights=np.zeros((n,4),float)
 for comp in range(components):
  for cls in [0,1,2]:
   vids=np.flatnonzero((labels==comp)&(klass==cls));nids=np.flatnonzero((nodeComponent==comp)&((nodeClass==cls)if cls<2 else np.ones(m,bool)))
   if not len(vids):continue
   assert len(nids)>0
   d,idx=cKDTree(nodes[nids]).query(p[vids],k=min(5,len(nids)));d=np.atleast_2d(d).reshape(len(vids),-1);idx=np.atleast_2d(idx).reshape(len(vids),-1);k=min(4,len(nids));dmax=d[:,-1]if len(nids)>4 else np.maximum(d[:,-1]*1.000001,h)
   weights=np.maximum(0,1-d[:,:k]/np.maximum(dmax[:,None],1e-20))**2;zero=weights.sum(1)==0;weights[zero,0]=1;weights/=weights.sum(1,keepdims=True);embedIds[vids,:k]=nids[idx[:,:k]];embedWeights[vids,:k]=weights
 G={'method':'source-embedded-arap-v1','nodes':nodes.ravel().tolist(),'nodeWeights':nodeW.tolist(),'nodeClass':nodeClass.tolist(),'edges':E.ravel().tolist(),'edgeWeights':ew.tolist(),'free':free.tolist(),'fixed':fixed.tolist(),'factor':{'lowerOffsets':lower.indptr.tolist(),'lowerIndices':lower.indices.tolist(),'lowerValues':lower.data.tolist(),'upperOffsets':upper.indptr.tolist(),'upperIndices':upper.indices.tolist(),'upperValues':upper.data.tolist(),'diagonal':diag.tolist()},'cellWorld':float(h),'unanchoredComponents':unanchored,'embeddingMethod':'4 nearest within source component; fixed vertices use their bone class; normalized compact weights'}
 (out/(id+'-graph.json')).write_text(json.dumps(G,separators=(',',':')));np.c_[embedIds.astype('<f4'),embedWeights.astype('<f4')].astype('<f4').tofile(out/(id+'-embedding.f32'));B['graphFile']=id+'-graph.json';B['graphEmbeddingFile']=id+'-embedding.f32';bindings.append(B);r={'id':id,'nodes':m,'freeNodes':len(free),'sourceVertices':n,'edges':len(E),'factorNonzeros':len(lower.data)+len(upper.data),'unanchored':unanchored};records.append(r);print(json.dumps(r),flush=True)
(out/'bindings.json').write_text(json.dumps(bindings,separators=(',',':')));(out/'report.json').write_text(json.dumps({'method':__doc__,'meshes':records},indent=2)+'\n')
