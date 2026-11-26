"""Exact rest-tangent response modes for the normalized source solid.
Basis construction is driven by bone-boundary affine motions, not muscle-name
coefficients. Corrections vanish at every original fixed attachment vertex.
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator,cg


class RestElasticity:
    def __init__(self,solid,bulk):
        self.solid=solid;self.bulk=bulk
        self.inverseT=solid.inverse.transpose(0,2,1)
        n=len(solid.free)*3
        def multiply(value):
            full=np.zeros_like(solid.rest);full[solid.free]=value.reshape(-1,3)
            return self.apply(full)[solid.free].ravel()
        self.operator=LinearOperator((n,n),matvec=multiply,dtype=float)
        self.preconditioner=LinearOperator((n,n),matvec=lambda value:solid.factor.solve(value.reshape(-1,3)).ravel(),dtype=float)

    def apply(self,displacement):
        s=self.solid;F=s.edges(displacement)@s.inverse
        trace=np.trace(F,axis1=1,axis2=2)
        stress=F+F.transpose(0,2,1)+(self.bulk-2/3)*trace[:,None,None]*np.eye(3)
        G=(stress@self.inverseT)*s.weights[:,None,None]
        out=np.zeros_like(displacement);np.add.at(out,s.tets[:,0],-G.sum(2))
        for k in range(3):np.add.at(out,s.tets[:,k+1],G[:,:,k])
        return out

    def solve_boundary(self,boundary,rtol=1e-7,maxiter=2000):
        s=self.solid;x=np.zeros_like(s.rest);x[s.fixed]=boundary
        rhs=-self.apply(x)[s.free].ravel();count=[0]
        def progress(_):count[0]+=1
        solution,info=cg(self.operator,rhs,M=self.preconditioner,rtol=rtol,atol=1e-14,maxiter=maxiter,callback=progress)
        x[s.free]=solution.reshape(-1,3)
        gradient=self.apply(x)[s.free]
        residual=float(np.linalg.norm(s.factor.solve(gradient),axis=1).max(initial=0))
        return x,{'cgInfo':int(info),'cgIterations':count[0],'preconditionedResidualMm':residual,'relativeLinearResidual':float(np.linalg.norm(gradient)/max(np.linalg.norm(rhs),1e-30))}


def build_modes(solid,lifting,pin_classes,bulk=30.,callback=None):
    solver=RestElasticity(solid,bulk);raw=[];reports=[]
    # Three translations and nine affine matrix entries span any rigid boundary
    # displacement R-I,t. They do not permit the fixed bone itself to stretch.
    for k in range(12):
        A=np.zeros((3,3));at=np.zeros(3)
        if k<3:at[k]=1
        else:A[(k-3)//3,(k-3)%3]=1
        boundary=np.zeros((len(solid.fixed),3));mask=pin_classes==1
        boundary[mask]=solid.rest[solid.fixed[mask]]@A.T+at
        response,report=solver.solve_boundary(boundary)
        harmonic=lifting[:,:3]@A.T+lifting[:,3,None]*at
        correction=response-harmonic
        if np.max(np.abs(correction[solid.fixed]))>1e-9:raise ValueError('Response mode changed a fixed source footprint')
        correction[solid.fixed]=0;raw.append(correction);reports.append({'mode':k,**report})
        if callback:callback(reports[-1])
    raw=np.stack(raw,axis=2)
    # Whiten in the coefficient-independent source FEM Dirichlet metric.
    gradients=np.stack([solid.edges(raw[:,:,k])@solid.inverse for k in range(raw.shape[2])],axis=3)
    weighted=(gradients*np.sqrt(solid.weights)[:,None,None,None]).reshape(-1,raw.shape[2]);S=weighted.T@weighted
    eigen,U=np.linalg.eigh(S);keep=eigen>eigen.max()*1e-12;transform=U[:,keep]/np.sqrt(eigen[keep]);basis=np.einsum('nik,km->nim',raw,transform)
    return {'vectorBasis':basis,'rawResponseModes':raw,'eigenvalues':eigen,'basisTransform':transform,'reports':reports,'discardedNumericalModes':int((~keep).sum())}
