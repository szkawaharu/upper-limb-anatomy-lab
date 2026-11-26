"""Matrix-free SPD strain/volume metric for the offline continuum reference.
This is an approximate Hessian, NOT the exact material Hessian. Armijo search
and exact element Jacobians still determine accepted steps. No muscle IDs.
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator,cg


def linearize(solid,x,bulk):
    F=solid.edges(x)@solid.inverse
    J=np.linalg.det(F)
    if np.min(J)<=0:raise ValueError('Non-positive input Jacobian')
    C=np.linalg.inv(F).transpose(0,2,1)
    scale=J**(-2/3)
    inverseT=solid.inverse.transpose(0,2,1)
    n=len(solid.free)*3
    def multiply(value):
        v=np.zeros_like(x);v[solid.free]=value.reshape(-1,3)
        dF=solid.edges(v)@solid.inverse
        divergence=np.einsum('tij,tij->t',C,dF)
        dP=scale[:,None,None]*dF+bulk*divergence[:,None,None]*C
        G=(dP@inverseT)*solid.weights[:,None,None]
        out=np.zeros_like(x);np.add.at(out,solid.tets[:,0],-G.sum(2))
        for k in range(3):np.add.at(out,solid.tets[:,k+1],G[:,:,k])
        return out[solid.free].ravel()
    return LinearOperator((n,n),matvec=multiply,dtype=float)


def minimize(solid,initial,bulk=3000.,iterations=30,cg_iterations=250,residual_tolerance_mm=.001,callback=None,preconditioner=None,preconditioner_info=None,cg_rtol=1e-3,preconditioner_factory=None):
    x=initial.copy();rows=[];reason='iteration-limit';n=len(solid.free)*3
    if preconditioner is None:
        preconditioner=LinearOperator((n,n),matvec=lambda v:solid.factor.solve(v.reshape(-1,3)).ravel(),dtype=float)
        preconditioner_info={'method':'source scalar FEM inverse, independent of bulk'}
    for iteration in range(iterations):
        energy,fullg,J=solid.evaluate(x,bulk);g=fullg[solid.free]
        residual=float(np.linalg.norm(solid.factor.solve(g),axis=1).max(initial=0))
        if residual<residual_tolerance_mm:reason='residual-tolerance';break
        if preconditioner_factory is not None:preconditioner=preconditioner_factory(x)
        H=linearize(solid,x,bulk);count=[0]
        def progress(_):count[0]+=1
        direction,info=cg(H,-g.ravel(),M=preconditioner,rtol=cg_rtol,atol=0,maxiter=cg_iterations,callback=progress)
        d=direction.reshape(-1,3);slope=float(np.sum(g*d))
        if not slope<0:reason='non-descent-cg-direction';break
        gain=1.;trial=x.copy();accepted=False
        for _ in range(40):
            trial[solid.free]=x[solid.free]+gain*d
            value,_,trialJ=solid.evaluate(trial,bulk,gradient=False)
            if np.isfinite(value)and value<=energy+.0001*gain*slope:accepted=True;break
            gain*=.5
        if not accepted:reason='line-search-stalled';break
        x=trial
        row={'iteration':iteration+1,'cgIterations':count[0],'cgInfo':int(info),'gain':gain,'energy':value,'volumeRatio':float(np.dot(solid.weights,trialJ)),'minJ':float(trialJ.min()),'maxJ':float(trialJ.max()),'previousResidualMaxMm':residual}
        rows.append(row)
        if callback:callback(row)
    energy,g,J=solid.evaluate(x,bulk);residual=float(np.linalg.norm(solid.factor.solve(g[solid.free]),axis=1).max(initial=0))
    if residual<residual_tolerance_mm:reason='residual-tolerance'
    return x,{'method':'matrix-free SPD strain/volume quasi-Newton; approximate Hessian','reason':reason,'converged':reason=='residual-tolerance','energy':energy,'bulkPenalty':bulk,'minJ':float(J.min()),'maxJ':float(J.max()),'JQuantiles':np.quantile(J,[0,.01,.5,.99,1]).tolist(),'volumeRatio':float(np.dot(solid.weights,J)),'nonPositiveElements':int((J<=0).sum()),'preconditionedGradientMaxMm':residual,'residualToleranceMm':residual_tolerance_mm,'linearPreconditioner':preconditioner_info,'cgRelativeTolerance':cg_rtol,'iterations':len(rows),'rows':rows,'activation':None,'force':None,'boneContactEvaluated':False}
