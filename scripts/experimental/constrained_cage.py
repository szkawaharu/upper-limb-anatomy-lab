"""Reference cage with recovered internal attachment facets.
Meshing is external; this locator contains no TetGen dependency. Exact boundary
facets are required by the generator rather than assuming unconstrained
Delaunay point interpolation preserves an anatomical attachment area.
"""
import numpy as np
from scipy.spatial import cKDTree
from rtree import index
from point_cage import PointCage

class ConstrainedCage(PointCage):
 def __init__(self,nodes,tets,classes):
  self.nodes=np.asarray(nodes);self.tets=np.asarray(tets);self.classes=np.asarray(classes)
  E=np.stack([self.nodes[self.tets[:,k]]-self.nodes[self.tets[:,0]]for k in [1,2,3]],axis=2)
  assert np.all(np.linalg.det(E)>0);self.inverse=np.linalg.inv(E);self.anchors=self.nodes[self.tets[:,0]];p=self.nodes[self.tets]
  self.index=index.Index(properties=index.Property(dimension=3))
  for i,(lo,hi)in enumerate(zip(p.min(1),p.max(1))):self.index.insert(i,(*lo,*hi))
  self.centroids=cKDTree(p.mean(1));self.node_tree=cKDTree(self.nodes);self.incident=np.full(len(nodes),-1,np.int32)
  for i,t in enumerate(self.tets):self.incident[t]=i
  self.neighbors=np.full((len(tets),4),-1,np.int32);faces={}
  for i,t in enumerate(self.tets):
   for opposite in range(4):
    key=tuple(sorted(np.delete(t,opposite)))
    if key in faces:
     j,k=faces.pop(key);self.neighbors[i,opposite]=j;self.neighbors[j,k]=i
    else:faces[key]=(i,opposite)

 def embed(self,points):
  points=np.asarray(points);_,ids=self.centroids.query(points);distance,nearest=self.node_tree.query(points);exact=distance==0
  if np.any(self.incident[nearest[exact]]<0):raise ValueError('Inactive exact source control')
  ids[exact]=self.incident[nearest[exact]];bary=np.zeros((len(points),4))
  for k in np.flatnonzero(exact):bary[k,np.flatnonzero(self.tets[ids[k]]==nearest[k])[0]]=1
  todo=np.flatnonzero(~exact);unresolved=[]
  for _ in range(64):
   if not len(todo):break
   tail=np.einsum('nij,nj->ni',self.inverse[ids[todo]],points[todo]-self.anchors[ids[todo]]);b=np.column_stack([1-tail.sum(1),tail]);bary[todo]=b
   bad=b.min(1)<-1e-9
   if not bad.any():todo=np.array([],int);break
   moving=todo[bad];face=np.argmin(b[bad],axis=1);next_ids=self.neighbors[ids[moving],face];valid=next_ids>=0;unresolved.extend(moving[~valid].tolist());ids[moving[valid]]=next_ids[valid];todo=moving[valid]
  unresolved=np.unique(np.r_[unresolved,todo]).astype(int)
  for k in unresolved:
   candidates=list(self.candidate_elements(points[k:k+1]));tail=np.einsum('nij,nj->ni',self.inverse[candidates],points[k]-self.anchors[candidates]);b=np.column_stack([1-tail.sum(1),tail]);score=b.min(1);best=int(np.argmax(score))
   if score[best]<-1e-9:raise ValueError('No valid constrained-cell embedding '+str(k)+': '+str(score[best]))
   ids[k]=candidates[best];bary[k]=b[best]
  for k in np.flatnonzero(exact):bary[k]=0;bary[k,np.flatnonzero(self.tets[ids[k]]==nearest[k])[0]]=1
  result=self.tets[ids]
  error=np.linalg.norm(self.positions(self.nodes,result,bary)-points,axis=1).max(initial=0)
  if error>1e-8:raise ValueError('Constrained embedding reconstruction error mm '+str(error))
  return result,bary,ids
