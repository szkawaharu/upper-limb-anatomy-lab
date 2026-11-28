"""Experimental element-wise inequality band on log J.
Numerical regularizer, NOT physical pressure or a calibrated muscle law.
All cells screened on every evaluation. Tiny cells receive a weight floor
only in the inequality term; the base material integration is unchanged."""
import numpy as np
from scipy.sparse import diags
from local_volume_constraints import log_volume_jacobian

class VolumeBand:
    def __init__(self,solid,rho,lower=.955,upper=1.045,multipliers=None):
        assert 0<lower<1<upper and np.isfinite(rho) and rho>0
        self.solid=solid;self.rho=float(rho)
        self.lower=np.log(lower);self.upper=np.log(upper)
        self.cells=np.arange(len(solid.tets))
        self.weights=np.maximum(solid.weights,1/len(self.cells))
        self.multipliers=np.zeros((2,len(self.cells))) if multipliers is None else np.array(multipliers,copy=True)
        assert self.multipliers.shape==(2,len(self.cells)) and np.isfinite(self.multipliers).all() and np.all(self.multipliers>=0)
        self.cached_nodes=None
    def state(self,x):
        if self.cached_nodes is None or not np.array_equal(self.cached_nodes,x):
            c,A=log_volume_jacobian(self.solid,x,self.cells)
            violations=np.stack([self.lower-c,c-self.upper])
            p=np.maximum(0,self.multipliers+self.rho*violations)
            self.cached_nodes=x.copy();self.cached=(c,A,p)
        return self.cached
    def evaluate(self,x,gradient=True):
        _,A,p=self.state(x)
        energy=float(np.sum(self.weights[None,:]*(p*p-self.multipliers*self.multipliers))/(2*self.rho))
        g=(A.T@(self.weights*(p[1]-p[0]))).reshape(x.shape) if gradient else None
        return energy,g
    def hessian(self,x):
        # Positive Gauss-Newton extra metric, not full constraint curvature.
        _,A,p=self.state(x)
        return A.T@diags(self.rho*self.weights*np.sum(p>0,axis=0))@A
    def update(self,x):
        return self.state(x)[2].copy()
