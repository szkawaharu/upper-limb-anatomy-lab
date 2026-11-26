"""Sparse source-rest strain/volume metric factorization, offline only.
This changes the linear solver, not the energy, boundaries or acceptance test.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu


def rest_metric_matrix(solid,bulk):
    gradients=np.concatenate([-solid.inverse.sum(1)[:,None,:],solid.inverse],axis=1)
    dot=np.einsum('tai,tbi->tab',gradients,gradients)
    # H_(a,i),(b,j) = w[(grad a . grad b) delta_ij + k grad a_i grad b_j].
    local=dot[:,:,None,:,None]*np.eye(3)[None,None,:,None,:]+bulk*np.einsum('tai,tbj->taibj',gradients,gradients)
    local*=solid.weights[:,None,None,None,None]
    ids=solid.tets[:,:,None]*3+np.arange(3)[None,None,:]
    shape=local.shape
    rows=np.broadcast_to(ids[:,:,:,None,None],shape).ravel();cols=np.broadcast_to(ids[:,None,None,:,:],shape).ravel()
    H=coo_matrix((local.ravel(),(rows,cols)),shape=(len(solid.rest)*3,)*2).tocsr()
    free=(solid.free[:,None]*3+np.arange(3)[None,:]).ravel()
    return H[free][:,free].tocsc()


def factor_rest_metric(solid,bulk):
    H=rest_metric_matrix(solid,bulk)
    factor=splu(H,permc_spec='MMD_AT_PLUS_A',diag_pivot_thresh=0)
    return factor,{'method':'source-rest coupled strain/volume metric sparse LU','matrixShape':list(H.shape),'matrixNonzeros':int(H.nnz),'factorNonzeros':int(factor.L.nnz+factor.U.nnz),'bulkPenalty':bulk}
