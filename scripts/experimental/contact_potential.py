"""Unilateral sampled contact potential on an embedded material surface.
Plane construction and actual bone verification are external. Source overlap
allowances are per material point; this never pushes an entire muscle by a
muscle-specific cosmetic coefficient or infers activation/force.
"""
import numpy as np
from scipy.sparse import coo_matrix,diags,csc_matrix

class PlaneContact:
 def __init__(self,node_count,vertices,barycentric,normals,offsets,weights,penalty):
  vertices=np.asarray(vertices);b=np.asarray(barycentric);n=np.asarray(normals);self.offsets=np.asarray(offsets);self.weights=np.asarray(weights);self.penalty=penalty;self.shape=(node_count,3)
  assert b.shape==vertices.shape and n.shape==(len(vertices),3)and self.weights.shape==(len(vertices),)
  assert np.all(self.weights>0)and penalty>0 and np.max(np.abs(np.linalg.norm(n,axis=1)-1),initial=0)<1e-8
  columns=(vertices[:,:,None]*3+np.arange(3)).reshape(len(vertices),12);values=(b[:,:,None]*n[:,None,:]).reshape(len(vertices),12);rows=np.broadcast_to(np.arange(len(vertices))[:,None],columns.shape)
  self.A=coo_matrix((values.ravel(),(rows.ravel(),columns.ravel())),shape=(len(vertices),node_count*3)).tocsr()

 def evaluate(self,x,gradient=True):
  violation=np.maximum(self.offsets-self.A@x.ravel(),0);energy=.5*self.penalty*np.dot(self.weights,violation**2)
  g=-(self.A.T@(self.penalty*self.weights*violation)).reshape(self.shape)if gradient else None
  return float(energy),g

 def hessian(self,x):
  active=self.offsets-self.A@x.ravel()>0
  return self.A.T@diags(self.penalty*self.weights*active)@self.A

class AugmentedSolid:
 def __init__(self,base,potentials):self.base=base;self.potentials=potentials
 def __getattr__(self,name):return getattr(self.base,name)
 def evaluate(self,x,bulk=30.,gradient=True):
  energy,g,J=self.base.evaluate(x,bulk,gradient)
  if not np.isfinite(energy):return energy,g,J
  for potential in self.potentials:
   value,d=potential.evaluate(x,gradient);energy+=value
   if gradient:g+=d
  return energy,g,J
 def extra_hessian(self,x):
  result=csc_matrix((x.size,x.size))
  for potential in self.potentials:result=result+potential.hessian(x)
  return result
