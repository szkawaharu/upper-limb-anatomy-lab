"""Offline vector-response reduced solid. Same exact element energy as the
full reference. Reduced convergence does not imply full FEM equilibrium, good
local volumes, collision-free anatomy or real-time performance.
"""
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from volumetric_reference import Solid


class VectorReducedSolid:
    def __init__(self,reference,tets,fixed,basis):
        self.solid=Solid(reference,tets,fixed);self.basis=np.asarray(basis)
        if self.basis.shape[:2] != self.solid.rest.shape:raise ValueError('Vector basis shape mismatch')
        if np.any(self.basis[fixed]!=0):raise ValueError('Response mode moves a fixed source attachment')
        self.G=np.stack([self.solid.edges(self.basis[:,:,k])@self.solid.inverse for k in range(self.basis.shape[2])],axis=3)
        weighted=(self.G*np.sqrt(self.solid.weights)[:,None,None,None]).reshape(-1,self.basis.shape[2])
        self.metric=weighted.T@weighted;self.metricFactor=cho_factor(self.metric)
        self.metricError=float(np.linalg.norm(self.metric-np.eye(len(self.metric))))

    def positions(self,base,q):return base+np.einsum('nik,k->ni',self.basis,q)

    def evaluate(self,base,q,bulk=30.,hessian=False):
        s=self.solid;F=s.edges(base)@s.inverse+np.einsum('tijk,k->tij',self.G,q)
        J=np.linalg.det(F)
        if not np.all(np.isfinite(J))or J.min()<=1e-9:return np.inf,None,J,None
        Ic=np.einsum('tij,tij->t',F,F);logJ=np.log(J);scale=J**(-2/3);C=np.linalg.inv(F).transpose(0,2,1)
        energy=float(np.dot(s.weights,.5*(scale*Ic-3)+.5*bulk*logJ**2))
        P=scale[:,None,None]*(F-Ic[:,None,None]*C/3)+bulk*logJ[:,None,None]*C
        gradient=np.einsum('t,tij,tijk->k',s.weights,P,self.G,optimize=True)
        H=None
        if hessian:
            weightedG=(self.G*np.sqrt(s.weights*scale)[:,None,None,None]).reshape(-1,len(q))
            divergence=np.einsum('tij,tijk->tk',C,self.G,optimize=True)*np.sqrt(s.weights)[:,None]
            H=weightedG.T@weightedG+bulk*(divergence.T@divergence)
        return energy,gradient,J,H

    def residual(self,g):
        displacement=np.einsum('nik,k->ni',self.basis,cho_solve(self.metricFactor,g))
        return float(np.linalg.norm(displacement,axis=1).max(initial=0))

    def minimize(self,base,q=None,bulk=30.,iterations=40,tolerance_mm=.001):
        q=np.zeros(self.basis.shape[2])if q is None else q.copy();rows=[];reason='iteration-limit'
        for iteration in range(iterations):
            energy,g,J,H=self.evaluate(base,q,bulk,hessian=True)
            if not np.isfinite(energy):raise ValueError('Reduced input is not positive-J feasible')
            residual=self.residual(g)
            if residual<tolerance_mm:reason='reduced-residual-tolerance';break
            d=cho_solve(cho_factor(H),-g);slope=float(np.dot(g,d))
            if not slope<0:reason='non-descent-reduced-direction';break
            gain=1.;accepted=False
            for _ in range(40):
                trial=q+gain*d;value,_,trialJ,_=self.evaluate(base,trial,bulk)
                if np.isfinite(value)and value<=energy+.0001*gain*slope:accepted=True;break
                gain*=.5
            if not accepted:reason='reduced-line-search-stalled';break
            q=trial;rows.append({'iteration':iteration+1,'gain':gain,'energy':value,'minJ':float(trialJ.min()),'maxJ':float(trialJ.max()),'previousProjectedResidualMm':residual})
        x=self.positions(base,q);energy,g,J,_=self.evaluate(base,q,bulk);_,fullg,_=self.solid.evaluate(x,bulk)
        projected=self.residual(g);full=float(np.linalg.norm(self.solid.factor.solve(fullg[self.solid.free]),axis=1).max(initial=0))
        if projected<tolerance_mm:reason='reduced-residual-tolerance'
        return x,q,{'method':'boundary-response vector reduced source solid','convergenceSpace':'reduced','converged':reason=='reduced-residual-tolerance','reason':reason,'iterations':len(rows),'energy':energy,'bulkPenalty':bulk,'volumeRatio':float(np.dot(self.solid.weights,J)),'minJ':float(J.min()),'maxJ':float(J.max()),'JQuantiles':np.quantile(J,[0,.01,.5,.99,1]).tolist(),'nonPositiveElements':int((J<=0).sum()),'projectedResidualMaxMm':projected,'fullResidualMaxMm':full,'fullReferenceConverged':full<tolerance_mm,'modes':len(q),'metricWhiteningError':self.metricError,'rows':rows,'force':None,'activation':None,'boneContactEvaluated':False}
