"""Common piecewise-affine embedding cage; preserves all source surface data.
A numerical background grid is not an anatomical solid. Attachment supports
are rigid Dirichlet patches; material quadrature is supplied independently.
"""
import itertools
import numpy as np

PERMUTATIONS=list(itertools.permutations(range(3)))

class GridCage:
 def __init__(self,points,spacing_mm):
  self.h=float(spacing_mm);self.origin=np.floor(np.min(points,axis=0)/self.h)*self.h-self.h
  self.cells=np.ceil((np.max(points,axis=0)-self.origin)/self.h).astype(int)+1;self.shape=self.cells+1
  grid=np.indices(self.shape).reshape(3,-1).T;self.nodes=self.origin+grid*self.h
  ids=np.indices(self.cells).reshape(3,-1).T;tets=[]
  for cell in ids:
   for perm in PERMUTATIONS:
    offsets=np.zeros((4,3),int)
    for k in range(3):offsets[k+1]=offsets[k];offsets[k+1,perm[k]]+=1
    tet=np.ravel_multi_index((cell+offsets).T,self.shape)
    p=self.nodes[tet]
    if np.linalg.det(np.column_stack([p[k]-p[0]for k in [1,2,3]]))<0:tet[[1,2]]=tet[[2,1]]
    tets.append(tet)
  self.tets=np.array(tets,np.int32)

 def embed(self,points):
  relative=(np.asarray(points)-self.origin)/self.h;cell=np.floor(relative).astype(int);fraction=relative-cell
  if np.any(cell<0)or np.any(cell>=self.cells):raise ValueError('Source points outside grid cage')
  permutation=np.argsort(-fraction,axis=1,kind='stable');order=np.array([PERMUTATIONS.index(tuple(row))for row in permutation]);tet_id=np.ravel_multi_index(cell.T,self.cells)*6+order
  vertices=self.tets[tet_id];p=self.nodes[vertices];E=np.stack([p[:,k]-p[:,0]for k in [1,2,3]],axis=2)
  tail=np.linalg.solve(E,(np.asarray(points)-p[:,0])[...,None])[...,0];bary=np.column_stack([1-tail.sum(1),tail]);reconstructed=np.einsum('ni,nij->nj',bary,p)
  if np.max(np.abs(reconstructed-points))>1e-9 or bary.min()<-1e-10:raise ValueError('Invalid source embedding')
  return vertices,bary,tet_id

 def fixed_patches(self,source_points,source_weights):
  vertices,bary,tet_id=self.embed(source_points);classes=np.full(len(self.nodes),-1,np.int8)
  for label in [0,1]:
   use=np.asarray(source_weights)==label
   selected=np.unique(vertices[use][bary[use]>1e-12])
   if np.any((classes[selected]>=0)&(classes[selected]!=label)):raise ValueError('Attachment supports overlap; refine the cage')
   classes[selected]=label
  return classes,vertices,bary,tet_id

 def material_weights(self,source_nodes,source_tets):
  # Positive four-point degree-two tetrahedron quadrature. The discontinuous
  # cage-cell density is approximated; source boundary volume is audited too.
  a=.5854101966249685;b=.1381966011250105;quadrature=np.full((4,4),b);np.fill_diagonal(quadrature,a)
  p=source_nodes[source_tets];E=np.stack([p[:,k]-p[:,0]for k in [1,2,3]],axis=2);volume=np.linalg.det(E)/6
  if np.any(volume<=0):raise ValueError('Material reference must have positive source elements')
  samples=np.einsum('qi,tij->tqj',quadrature,p).reshape(-1,3);_,_,ids=self.embed(samples)
  occupied=np.bincount(ids,weights=np.repeat(volume/4,4),minlength=len(self.tets));physical=occupied/occupied.sum()
  # Empty background has no assigned muscle mass. A tiny source-independent
  # numerical regularizer removes null modes; this is separately disclosed.
  regularizer=float(physical.max()*1e-6);weights=physical+regularizer
  return weights/weights.sum(),{'emptyElementRegularizer':regularizer,'occupiedElements':int((physical>0).sum()),'backgroundElements':len(self.tets),'sourceVolumeMm3':float(volume.sum()),'quadrature':'positive symmetric 4-point source tetrahedra; cell-density approximation'}

 def positions(self,nodes,vertices,bary):return np.einsum('ni,nij->nj',bary,nodes[vertices])
