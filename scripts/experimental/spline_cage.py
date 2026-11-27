"""C1 quadratic Cartesian B-spline embedding. No muscle-local curve/frame.
All source vertices and faces are retained. This is a geometric field, not
measured muscle activity or measured material mechanics.
"""
import itertools
import numpy as np
from scipy.sparse import coo_matrix,csr_matrix
from scipy.linalg import cho_factor,cho_solve
from scipy.sparse.linalg import splu
from spline_regularization import BendingExtension
from constitutive_metric import tangent
OFFSETS=np.array(list(itertools.product([-1,0,1],repeat=3)))

class SplineCage:
 def __init__(self,points,spacing_mm):
  self.h=float(spacing_mm);self.origin=np.floor(np.min(points,axis=0)/self.h)*self.h-2*self.h
  self.shape=np.ceil((np.max(points,axis=0)-self.origin)/self.h).astype(int)+3
  self.nodes=self.origin+np.indices(self.shape).reshape(3,-1).T*self.h

 def embed(self,points):
  r=(np.asarray(points)-self.origin)/self.h;i=np.floor(r+.5).astype(int);u=r-i
  indices=i[:,None,:]+OFFSETS[None,:,:]
  if np.any(indices<0)or np.any(indices>=self.shape):raise ValueError('Spline support outside cage')
  ids=np.ravel_multi_index(indices.transpose(2,0,1),self.shape)
  values=np.stack([.5*(.5-u)**2,.75-u*u,.5*(.5+u)**2],axis=2)
  derivative=np.stack([u-.5,-2*u,u+.5],axis=2)/self.h
  w=[];g=[]
  for off in OFFSETS+1:
   v=[values[:,k,off[k]]for k in range(3)];d=[derivative[:,k,off[k]]for k in range(3)]
   w.append(v[0]*v[1]*v[2]);g.append(np.column_stack([d[0]*v[1]*v[2],v[0]*d[1]*v[2],v[0]*v[1]*d[2]]))
  return ids,np.array(w).T,np.array(g).transpose(1,0,2)

 def positions(self,nodes,embedding):
  ids,w,_=embedding;return np.einsum('qi,qij->qj',w,nodes[ids])

 def gradients(self,nodes,embedding):
  ids,_,g=embedding;return np.einsum('qai,qaj->qij',nodes[ids],g)

 def fixed_patches(self,points,weights):
  ids,_,_=self.embed(points);classes=np.full(len(self.nodes),-1,np.int8)
  for label in [0,1]:
   support=np.unique(ids[np.asarray(weights)==label])
   if np.any((classes[support]>=0)&(classes[support]!=label)):raise ValueError('Spline attachment supports overlap; refine spatial resolution')
   classes[support]=label
  return classes

 def quadrature(self,source_nodes,source_tets):
  a=.5854101966249685;b=.1381966011250105;Q=np.full((4,4),b);np.fill_diagonal(Q,a)
  p=source_nodes[source_tets];volume=np.einsum('ti,ti->t',p[:,1]-p[:,0],np.cross(p[:,2]-p[:,0],p[:,3]-p[:,0]))/6
  if np.any(volume<=0):raise ValueError('Non-positive source material element')
  samples=np.einsum('qi,tij->tqj',Q,p).reshape(-1,3);weights=np.repeat(volume/4,4)
  bins=np.floor((samples-self.origin)/(self.h/2)).astype(int);unique,inverse=np.unique(bins,axis=0,return_inverse=True);mass=np.bincount(inverse,weights=weights);centres=np.column_stack([np.bincount(inverse,weights=weights*samples[:,k])/mass for k in range(3)])
  return centres,mass/mass.sum(),{'sourceQuadraturePoints':len(samples),'integrationClusters':len(centres),'integrationCellMm':self.h/2,'sourceVolumeMm3':float(volume.sum()),'method':'positive 4-point source quadrature aggregated by half-control-spacing cells'}


class SplineSolid:
 def __init__(self,cage,points,weights,classes,regularization_ratio=1e-6):
  self.cage=cage;self.rest=cage.nodes;self.classes=classes;self.fixed=np.flatnonzero(classes>=0);self.integration=cage.embed(points);self.weights=np.asarray(weights)/np.sum(weights);self.ids,self.basis,self.gradient=self.integration
  active=np.unique(self.ids);self.free=np.setdiff1d(active,self.fixed);self.inactive=np.setdiff1d(np.arange(len(self.rest)),np.union1d(active,self.fixed))
  # Nodes influencing the quadrature define all free DOFs; no arbitrary empty
  # background material is introduced to regularize unused nodes.
  g=self.gradient;local=np.einsum('qai,qbi,q->qab',g,g,self.weights);rows=np.broadcast_to(self.ids[:,:,None],local.shape).ravel();cols=np.broadcast_to(self.ids[:,None,:],local.shape).ravel();L=coo_matrix((local.ravel(),(rows,cols)),shape=(len(self.rest),)*2).tocsr();self.bending=BendingExtension(cage,np.union1d(active,self.fixed));A=self.bending.active;self.K=self.bending.condensed;free_local=self.bending.index[self.free]
  primary=L[self.free][:,self.free].toarray();prior=self.K[np.ix_(free_local,free_local)];scale=np.linalg.eigvalsh(primary).max()/np.linalg.eigvalsh(prior).max();self.bending_weight=regularization_ratio*scale
  correction=coo_matrix((self.K.ravel(),(np.repeat(A,len(A)),np.tile(A,len(A)))),shape=L.shape).tocsr();L=L+self.bending_weight*correction
  self.factor=splu(L[self.free][:,self.free].tocsc());self.coupling=L[self.free][:,self.fixed]
  self.free_gradients=np.zeros((len(self.ids),3,len(self.free)));free_index=np.full(len(self.rest),-1,int);free_index[self.free]=np.arange(len(self.free))
  for j in range(self.ids.shape[1]):
   selected=free_index[self.ids[:,j]];use=np.flatnonzero(selected>=0);self.free_gradients[use,:,selected[use]]=self.gradient[use,j,:]
  self.free_bending=np.kron(self.bending_weight*self.K[np.ix_(free_local,free_local)],np.eye(3))
  self.regularization={'relativeStiffness':regularization_ratio,'bendingWeight':self.bending_weight,'primaryEigenvalueRange':np.linalg.eigvalsh(primary)[[0,-1]].tolist(),'regularizedEigenvalueRange':np.linalg.eigvalsh(L[self.free][:,self.free].toarray())[[0,-1]].tolist(),'interpretation':'numerical minimum-bending prior, no anatomical mass'}

 def predict(self,x,target):
  out=x.copy();out[self.free]-=self.factor.solve(self.coupling@(target-x[self.fixed]));out[self.fixed]=target;return self.bending.extend(out,self.rest)

 def evaluate(self,x,bulk=30.,hessian=False):
  F=self.cage.gradients(x,self.integration);J=np.linalg.det(F)
  if not np.all(np.isfinite(J))or J.min()<=1e-9:return np.inf,None,J,None
  if hasattr(self,'guards'):
   guard_J=np.linalg.det(self.cage.gradients(x,self.guards))
   if not np.all(np.isfinite(guard_J))or guard_J.min()<=1e-9:return np.inf,None,J,None
  Ic=np.einsum('qij,qij->q',F,F);logJ=np.log(J);scale=J**(-2/3);C=np.linalg.inv(F).transpose(0,2,1)
  u=x[self.bending.active]-self.rest[self.bending.active];bend_g=self.bending_weight*(self.K@u);bend_energy=.5*float(np.sum(u*bend_g))
  E=bend_energy+float(np.dot(self.weights,.5*(scale*Ic-3)+.5*bulk*logJ**2));P=scale[:,None,None]*(F-Ic[:,None,None]*C/3)+bulk*logJ[:,None,None]*C
  local=np.einsum('qij,qaj,q->qai',P,self.gradient,self.weights);g=np.zeros_like(x);np.add.at(g,self.ids,local);g[self.bending.active]+=bend_g
  H=None
  if hessian:
   root,self.metric_info=tangent(F,bulk,project=True)
   # D_(i,j),(node,k) = delta_(i,k) gradN_(j,node).
   D=np.einsum('ik,qjn->qijnk',np.eye(3),self.free_gradients).reshape(len(J),9,-1)
   W=np.einsum('qir,qik->qrk',root,D,optimize=True)*np.sqrt(self.weights)[:,None,None]
   W=W.reshape(-1,len(self.free)*3);H=W.T@W+self.free_bending
  return E,g,J,H

 def minimize(self,initial,bulk=30.,iterations=100,tolerance_mm=.001,callback=None):
  x=self.bending.extend(initial,self.rest);rows=[];reason='iteration-limit'
  for k in range(iterations):
   E,g,J,H=self.evaluate(x,bulk,hessian=True)
   if not np.isfinite(E):raise ValueError('Non-positive spline integration Jacobian')
   residual=float(np.linalg.norm(self.factor.solve(g[self.free]),axis=1).max(initial=0))
   if residual<tolerance_mm:reason='residual-tolerance';break
   factor=cho_factor(H);d=cho_solve(factor,-g[self.free].ravel()).reshape(-1,3);slope=float(np.sum(g[self.free]*d));gain=1.;success=False;trial=x.copy()
   if not slope<0:reason='non-descent';break
   for _ in range(40):
    trial[self.free]=x[self.free]+gain*d;trial=self.bending.extend(trial,self.rest);value,_,nextJ,_=self.evaluate(trial,bulk)
    if np.isfinite(value)and value<=E+.0001*gain*slope:success=True;break
    gain*=.5
   if not success:reason='line-search-stalled';break
   x=trial;rows.append({'iteration':k+1,'energy':value,'gain':gain,'minJ':float(nextJ.min()),'maxJ':float(nextJ.max()),'previousResidualMm':residual})
   if callback:callback(rows[-1])
  E,g,J,_=self.evaluate(x,bulk);residual=float(np.linalg.norm(self.factor.solve(g[self.free]),axis=1).max(initial=0))
  if residual<tolerance_mm:reason='residual-tolerance'
  return x,{'converged':reason=='residual-tolerance','reason':reason,'iterations':len(rows),'residualMaxMm':residual,'energy':E,'JQuantiles':np.quantile(J,[0,.01,.5,.99,1]).tolist(),'nonPositiveElements':int((J<=0).sum()),'quadratureVolumeRatio':float(np.dot(self.weights,J)),'rows':rows,'regularization':self.regularization,'optimizationMetric':getattr(self,'metric_info',None)}
