"""Sparse source-rest strain/volume metric factorization, offline only.
This changes the linear solver, not the energy, boundaries or acceptance test.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu,LinearOperator


def rest_metric_matrix(solid,bulk):
    return strain_metric_matrix(solid,solid.rest,bulk)


def strain_metric_matrix(solid,x,bulk):
    """Assemble the same current SPD metric as the independent matrix-free
    operator. It is an approximate material Hessian, not an energy change.
    """
    gradients=np.concatenate([-solid.inverse.sum(1)[:,None,:],solid.inverse],axis=1)
    F=solid.edges(x)@solid.inverse;J=np.linalg.det(F)
    if J.min()<=0:raise ValueError('Non-positive metric input')
    C=np.linalg.inv(F).transpose(0,2,1)
    volume_grad=np.einsum('tij,taj->tai',C,gradients)
    dot=np.einsum('tai,tbi->tab',gradients,gradients)
    # H_(a,i),(b,j) = w[J^-2/3(grad a . grad b) delta_ij
    #                     + k (F^-T grad a)_i (F^-T grad b)_j].
    local=(dot*J[:,None,None]**(-2/3))[:,:,None,:,None]*np.eye(3)[None,None,:,None,:]+bulk*np.einsum('tai,tbj->taibj',volume_grad,volume_grad)
    local*=solid.weights[:,None,None,None,None]
    ids=solid.tets[:,:,None]*3+np.arange(3)[None,None,:]
    shape=local.shape
    rows=np.broadcast_to(ids[:,:,:,None,None],shape).ravel();cols=np.broadcast_to(ids[:,None,None,:,:],shape).ravel()
    H=coo_matrix((local.ravel(),(rows,cols)),shape=(len(solid.rest)*3,)*2).tocsr()
    free=(solid.free[:,None]*3+np.arange(3)[None,:]).ravel()
    return H[free][:,free].tocsc()


def factor_rest_metric(solid,bulk):
    return factor_strain_metric(solid,solid.rest,bulk)


def factor_strain_metric(solid,x,bulk):
    H=strain_metric_matrix(solid,x,bulk)
    factor=splu(H,permc_spec='MMD_AT_PLUS_A',diag_pivot_thresh=0)
    return factor,{'method':'frozen-pose coupled strain/volume metric sparse LU','matrixShape':list(H.shape),'matrixNonzeros':int(H.nnz),'factorNonzeros':int(factor.L.nnz+factor.U.nnz),'bulkPenalty':bulk}


def warped_rest_preconditioner(solid,x,factor):
    """Transport the linear preconditioner by source-volume-averaged nodal
    polar rotations. This is a congruence, so symmetry/positive definiteness
    survive. It does not rotate the mesh or alter the nonlinear strain energy.
    A common rigid motion transports the full rest metric exactly.
    """
    F=solid.edges(x)@solid.inverse
    U,_,Vt=np.linalg.svd(F);sign=np.linalg.det(U@Vt)
    U[:,:,2]*=sign[:,None];rotations=U@Vt
    sums=np.zeros((len(x),3,3));mass=np.zeros(len(x))
    for k in range(4):
        np.add.at(sums,solid.tets[:,k],rotations*solid.weights[:,None,None])
        np.add.at(mass,solid.tets[:,k],solid.weights)
    active=np.flatnonzero(mass);averages=sums[active]/mass[active,None,None]
    U,_,Vt=np.linalg.svd(averages);U[:,:,2]*=np.linalg.det(U@Vt)[:,None]
    nodal=np.broadcast_to(np.eye(3),(len(x),3,3)).copy();nodal[active]=U@Vt;R=nodal[solid.free]
    def inverse(value):
        local=np.einsum('nji,nj->ni',R,value.reshape(-1,3))
        solved=factor.solve(local.ravel()).reshape(-1,3)
        return np.einsum('nij,nj->ni',R,solved).ravel()
    n=len(solid.free)*3
    return LinearOperator((n,n),matvec=inverse,dtype=float)
