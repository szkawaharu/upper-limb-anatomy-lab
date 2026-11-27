"""Sparse projected material tangent for the piecewise-affine source cage.
Local eigenvalue projection changes the optimization metric only.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from constitutive_metric import tangent

def assemble(solid,x,bulk,project=True):
 F=solid.edges(x)@solid.inverse
 gradients=np.concatenate([-solid.inverse.sum(1)[:,None,:],solid.inverse],axis=1)
 D=np.einsum('ik,taj->tijak',np.eye(3),gradients).reshape(len(F),9,12)
 if project:
  root,info=tangent(F,bulk,True);W=np.einsum('tir,tik->trk',root,D,optimize=True)*np.sqrt(solid.weights)[:,None,None];H=np.einsum('tri,trj->tij',W,W)
 else:
  local=tangent(F,bulk,False);H=np.einsum('tir,tij,tjs->trs',D,local,D,optimize=True)*solid.weights[:,None,None];info={'projection':'none; exact assembled global tangent, direction must be descent'}
 ids=(solid.tets[:,:,None]*3+np.arange(3)).reshape(-1,12);rows=np.broadcast_to(ids[:,:,None],H.shape).ravel();cols=np.broadcast_to(ids[:,None,:],H.shape).ravel();matrix=coo_matrix((H.ravel(),(rows,cols)),shape=(len(x)*3,)*2).tocsr()
 if hasattr(solid,'extra_hessian'):matrix=matrix+solid.extra_hessian(x)
 free=(solid.free[:,None]*3+np.arange(3)).ravel();matrix=matrix[free][:,free].tocsc()
 return matrix,info

def factor(solid,x,bulk,project=True):
 matrix,info=assemble(solid,x,bulk,project)
 linear=splu(matrix,permc_spec='MMD_AT_PLUS_A',diag_pivot_thresh=0,options={'SymmetricMode':True})
 matrix_norm=float(np.asarray(abs(matrix).sum(axis=1)).max())
 info['factorization']='symmetric ordering, zero diagonal pivot threshold; checked solve, descent and line-search safeguards'
 class CheckedFactor:
  def solve(self,rhs):
   step=linear.solve(rhs);denominator=matrix_norm*np.max(np.abs(step))+np.max(np.abs(rhs));error=float(np.max(np.abs(matrix@step-rhs))/max(denominator,1e-300));info['linearBackwardError']=error
   if not np.isfinite(step).all()or not np.isfinite(error)or error>1e-10:raise RuntimeError('Linear tangent solve failed backward-error check')
   return step
 return CheckedFactor(),info
