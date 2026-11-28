"""Offline principal-stretch regularizer; not a physiological strain limit.
Penalize log singular values outside a diagnostic band. Exact scalar energy
and gradient, PSD Gauss-Newton search metric. Base material weights unchanged.
Repeated singular values have identical scalar derivatives; the GN metric
is only a search metric, not the full spectral Hessian.
"""
import numpy as np
from scipy.sparse import coo_matrix,diags

class StretchBand:
    def __init__(self,solid,rho,lower=.5,upper=2.):
        assert 0<lower<1<upper and np.isfinite(rho) and rho>0
        self.solid=solid;self.rho=float(rho)
        self.lower=np.log(lower);self.upper=np.log(upper)
        self.weights=np.maximum(solid.weights,1/len(solid.tets))
        self.gradients=np.concatenate([-solid.inverse.sum(1)[:,None,:],solid.inverse],axis=1)
        self.cached_nodes=None
    def state(self,x):
        if self.cached_nodes is None or not np.array_equal(self.cached_nodes,x):
            F=self.solid.edges(x)@self.solid.inverse
            U,s,Vt=np.linalg.svd(F)
            assert np.isfinite(s).all() and s.min()>0
            c=np.log(s);e=c-np.clip(c,self.lower,self.upper)
            # d log(s_i)/d F = u_i v_i^T / s_i.
            active=np.nonzero(e)
            k,j=active
            deriv=np.einsum('ni,nj->nij',U[k,:,j],Vt[k,j,:])/s[k,j,None,None]
            local=np.einsum('nij,naj->nai',deriv,self.gradients[k])
            cols=(self.solid.tets[k,:,None]*3+np.arange(3)).reshape(-1,12)
            rows=np.broadcast_to(np.arange(len(k))[:,None],cols.shape)
            A=coo_matrix((local.ravel(),(rows.ravel(),cols.ravel())),shape=(len(k),x.size)).tocsr()
            self.cached_nodes=x.copy();self.cached=(s,e,A,k,j)
        return self.cached
    def evaluate(self,x,gradient=True):
        s,e,A,k,j=self.state(x)
        value=.5*self.rho*np.sum(self.weights[:,None]*e*e)
        g=(A.T@(self.rho*self.weights[k]*e[k,j])).reshape(x.shape) if gradient else None
        return float(value),g
    def hessian(self,x):
        _,_,A,k,_=self.state(x)
        return A.T@diags(self.rho*self.weights[k])@A
