"""Sparse projected material tangent for the piecewise-affine source cage.
Local eigenvalue projection changes the optimization metric only.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from constitutive_metric import tangent

def assemble(solid,x,bulk):
 F=solid.edges(x)@solid.inverse;root,info=tangent(F,bulk,True)
 gradients=np.concatenate([-solid.inverse.sum(1)[:,None,:],solid.inverse],axis=1)
 D=np.einsum('ik,taj->tijak',np.eye(3),gradients).reshape(len(F),9,12)
 W=np.einsum('tir,tik->trk',root,D,optimize=True)*np.sqrt(solid.weights)[:,None,None]
 H=np.einsum('tri,trj->tij',W,W);ids=(solid.tets[:,:,None]*3+np.arange(3)).reshape(-1,12);rows=np.broadcast_to(ids[:,:,None],H.shape).ravel();cols=np.broadcast_to(ids[:,None,:],H.shape).ravel();matrix=coo_matrix((H.ravel(),(rows,cols)),shape=(len(x)*3,)*2).tocsr();free=(solid.free[:,None]*3+np.arange(3)).ravel();matrix=matrix[free][:,free].tocsc()
 return matrix,info

def factor(solid,x,bulk):
 matrix,info=assemble(solid,x,bulk)
 return splu(matrix,permc_spec='MMD_AT_PLUS_A',diag_pivot_thresh=0),info
