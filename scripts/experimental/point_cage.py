"""A neutral unstructured control cage with exact attachment vertices.
Grid support pinning enlarges a Dirichlet patch into adjacent tissue. Here every
declared attachment point is a control vertex; no other control is fixed on its
behalf. The original anatomy is embedded, not replaced by the cage hull.
"""
import numpy as np
from scipy.spatial import Delaunay,cKDTree
from rtree import index

class PointCage:
 def __init__(self,points,weights,spacing):
  self.source=np.asarray(points,float);w=np.asarray(weights);self.fixed_source=np.flatnonzero(np.isin(w,[0,1]));self.h=spacing
  low=np.floor(self.source.min(0)/spacing)*spacing-spacing;high=np.ceil(self.source.max(0)/spacing)*spacing+spacing
  grid=np.stack(np.meshgrid(*[np.arange(low[k],high[k]+spacing*.1,spacing)for k in range(3)],indexing='ij'),axis=-1).reshape(-1,3)
  attachment=self.source[self.fixed_source];distance=cKDTree(attachment).query(grid)[0]
  # Only exact duplicate grid points are omitted; no tissue or footprint
  # merging, user-tuned pin radius, joggle, or per-muscle coefficient.
  grid=grid[distance>np.finfo(float).eps*max(1.,np.abs(grid).max())*16]
  self.nodes=np.concatenate([attachment,grid]);self.classes=np.r_[w[self.fixed_source],np.full(len(grid),-1)].astype(np.int8)
  self.triangulation=Delaunay(self.nodes);self.tets=self.triangulation.simplices.copy()
  assert np.array_equal(np.unique(self.tets),np.arange(len(self.nodes))),'A control/attachment point was omitted'
  E=np.stack([self.nodes[self.tets[:,k]]-self.nodes[self.tets[:,0]]for k in [1,2,3]],axis=2);det=np.linalg.det(E)
  # Qhull can emit exactly coplanar lattice simplices. They enclose no
  # volume; remove only exact zeros, without moving source points or dropping
  # any positive-volume sliver to improve a deformation score.
  self.removed_zero_simplices=int((det==0).sum());keep=det!=0;self.simplex_map=np.full(len(det),-1,np.int32);self.simplex_map[keep]=np.arange(keep.sum());self.tets=self.tets[keep];det=det[keep]
  assert np.isin(np.arange(len(attachment)),np.unique(self.tets)).all(),'An attachment became inactive'
  reverse=det<0;self.tets[reverse,1],self.tets[reverse,2]=self.tets[reverse,2].copy(),self.tets[reverse,1].copy()
  self.inverse=np.linalg.inv(np.stack([self.nodes[self.tets[:,k]]-self.nodes[self.tets[:,0]]for k in [1,2,3]],axis=2));self.anchors=self.nodes[self.tets[:,0]]
  self.index=index.Index(properties=index.Property(dimension=3))
  vertices=self.nodes[self.tets]
  for i,(lo,hi)in enumerate(zip(vertices.min(1),vertices.max(1))):self.index.insert(i,(*lo,*hi))

 def candidate_elements(self,triangle):
  lo=triangle.min(0);hi=triangle.max(0);eps=np.finfo(float).eps*max(1.,np.abs(triangle).max())*32
  return self.index.intersection((* (lo-eps),* (hi+eps)))

 def embed(self,points):
  points=np.asarray(points);original=self.triangulation.find_simplex(points,tol=1e-9)
  if np.any(original<0):raise ValueError('Source point outside the control cage')
  ids=self.simplex_map[original]
  invalid=ids<0;ids[invalid]=0
  tail=np.einsum('nij,nj->ni',self.inverse[ids],points-self.anchors[ids]);bary=np.column_stack([1-tail.sum(1),tail])
  # The walking locator can pick a numerically ill-conditioned sliver at a
  # boundary. Search geometrically overlapping positive cells instead of
  # accepting a negative weight or moving an anatomical source point.
  for k in np.flatnonzero(invalid|(bary.min(1)<-1e-9)):
   candidates=list(self.candidate_elements(points[k:k+1]));tail=np.einsum('nij,nj->ni',self.inverse[candidates],points[k]-self.anchors[candidates]);b=np.column_stack([1-tail.sum(1),tail]);score=b.min(1);best=int(np.argmax(score))
   if score[best]<-1e-9:raise ValueError('No positive-volume embedding for point '+str(k)+'; best barycentric minimum '+str(score[best]))
   ids[k]=candidates[best];bary[k]=b[best]
  # At declared fixed vertices, use the exact Kronecker coordinate, avoiding
  # conditioning loss in a nearly coplanar boundary control tetrahedron.
  dist,nearest=cKDTree(self.nodes).query(points);exact=dist==0
  for i in np.flatnonzero(exact):
   slots=np.flatnonzero(self.tets[ids[i]]==nearest[i])
   if len(slots):bary[i]=0;bary[i,slots[0]]=1
  return self.tets[ids],bary,ids

 def fixed_patches(self,points,weights):
  ids,bary,tets=self.embed(points)
  # Subdivision vertices on a fixed surface triangle are fixed by their
  # barycentric interpolation of the ORIGINAL fixed control vertices.
  for end in [0,1]:
   for i in np.flatnonzero(np.asarray(weights)==end):
    support=ids[i][bary[i]>1e-9]
    if np.any(self.classes[support]!=end):raise ValueError('Fixed surface vertex '+str(i)+' at end '+str(end)+' is not represented by its own fixed controls; support classes '+str(self.classes[support].tolist()))
  return self.classes.copy(),ids,bary,tets

 def positions(self,nodes,vertices,bary):return np.einsum('ni,nij->nj',bary,nodes[vertices])

 def material_weights(self,source_nodes,source_tets):
  q=np.full((4,4),.1381966011250105);np.fill_diagonal(q,.5854101966249685);p=source_nodes[source_tets];E=np.stack([p[:,k]-p[:,0]for k in [1,2,3]],axis=2);volume=np.linalg.det(E)/6
  assert np.all(volume>0);samples=np.einsum('qi,tij->tqj',q,p).reshape(-1,3);_,_,ids=self.embed(samples);physical=np.bincount(ids,weights=np.repeat(volume/4,4),minlength=len(self.tets));physical/=physical.sum();regularizer=float(physical.max()*1e-6);weights=physical+regularizer
  return weights/weights.sum(),{'emptyElementRegularizer':regularizer,'occupiedElements':int((physical>0).sum()),'backgroundElements':len(self.tets),'sourceVolumeMm3':float(volume.sum()),'quadrature':'positive symmetric 4-point source tetrahedra; cell-density approximation'}
