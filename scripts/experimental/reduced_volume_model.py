"""Reduced source solid prototype. Uses the same full-element energy and
positive-J acceptance as the continuum reference. Convergence in a reduced
space is reported separately from the full FEM residual. Not a web runtime.
"""
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from volumetric_reference import Solid

def affine_boundary_lift(reference,channels,A,at,B,bt):
    """Common rigid motion and both original fixed footprints are exact."""
    return reference@A.T+at+channels[:,:3]@(B-A).T+channels[:,3,None]*(bt-at)


class ReducedSolid:
    def __init__(self,reference,tets,fixed,basis):
        self.solid=Solid(reference,tets,fixed)
        self.basis=basis
        inverse=self.solid.inverse
        gradients=np.concatenate([-inverse.sum(1)[:,None,:],inverse],axis=1)
        self.G=np.einsum('tij,tik->tjk',gradients,basis[tets])
        weighted=self.G*np.sqrt(self.solid.weights)[:,None,None]
        gram=weighted.reshape(-1,basis.shape[1]).T@weighted.reshape(-1,basis.shape[1])
        self.metricError=float(np.linalg.norm(gram-np.eye(len(gram))))
        self.metricFactor=cho_factor(gram)
        if np.any(basis[fixed]!=0):raise ValueError('Reduced mode moves a fixed attachment')

    def evaluate(self,base,q,bulk=30.,hessian=False):
        x=base+self.basis@q
        energy,g,J=self.solid.evaluate(x,bulk)
        if not np.isfinite(energy):return energy,None,J,None
        gradient=self.basis.T@g
        H=None
        if hessian:
            F=self.solid.edges(x)@self.solid.inverse
            C=np.linalg.inv(F).transpose(0,2,1)
            B=np.einsum('tji,tim->tmj',C,self.G,optimize=True).reshape(len(F),-1)
            weightedB=B*np.sqrt(self.solid.weights)[:,None]
            weightedG=(self.G*np.sqrt(self.solid.weights*J**(-2/3))[:,None,None]).reshape(-1,self.basis.shape[1])
            H=np.kron(weightedG.T@weightedG,np.eye(3))+bulk*(weightedB.T@weightedB)
        return energy,gradient,J,H

    def minimize(self,base,q=None,bulk=30.,iterations=30,tolerance_mm=.001,callback=None):
        if q is None:q=np.zeros((self.basis.shape[1],3))
        else:q=q.copy()
        rows=[];reason='iteration-limit'
        for k in range(iterations):
            energy,g,J,H=self.evaluate(base,q,bulk,hessian=True)
            if not np.isfinite(energy):raise ValueError('Reduced initial pose is not positive-J feasible')
            residual=float(np.linalg.norm(self.basis@cho_solve(self.metricFactor,g),axis=1).max(initial=0))
            if residual<tolerance_mm:reason='reduced-residual-tolerance';break
            d=cho_solve(cho_factor(H),-g.ravel()).reshape(q.shape)
            slope=float(np.sum(g*d))
            if not slope<0:reason='non-descent-reduced-direction';break
            gain=1.;accepted=False
            for _ in range(40):
                trial=q+gain*d;value,_,newJ,_=self.evaluate(base,trial,bulk)
                if np.isfinite(value)and value<=energy+.0001*gain*slope:accepted=True;break
                gain*=.5
            if not accepted:reason='reduced-line-search-stalled';break
            q=trial;row={'iteration':k+1,'energy':value,'gain':gain,'minJ':float(newJ.min()),'maxJ':float(newJ.max()),'volumeRatio':float(np.dot(self.solid.weights,newJ)),'previousProjectedResidualMm':residual};rows.append(row)
            if callback:callback(row)
        x=base+self.basis@q;energy,g,J=self.solid.evaluate(x,bulk);projected=self.basis.T@g
        residual=float(np.linalg.norm(self.basis@cho_solve(self.metricFactor,projected),axis=1).max(initial=0))
        full_residual=float(np.linalg.norm(self.solid.factor.solve(g[self.solid.free]),axis=1).max(initial=0))
        if residual<tolerance_mm:reason='reduced-residual-tolerance'
        return x,q,{'method':'continuous Dirichlet affine-patch reduced source solid','convergenceSpace':'reduced','converged':reason=='reduced-residual-tolerance','reason':reason,'iterations':len(rows),'energy':energy,'bulkPenalty':bulk,'volumeRatio':float(np.dot(self.solid.weights,J)),'minJ':float(J.min()),'maxJ':float(J.max()),'JQuantiles':np.quantile(J,[0,.01,.5,.99,1]).tolist(),'nonPositiveElements':int((J<=0).sum()),'projectedResidualMaxMm':residual,'fullResidualMaxMm':full_residual,'fullReferenceConverged':full_residual<tolerance_mm,'modes':len(q),'metricWhiteningError':self.metricError,'rows':rows,'force':None,'activation':None,'boneContactEvaluated':False}
