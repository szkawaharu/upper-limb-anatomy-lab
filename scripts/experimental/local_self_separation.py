"""Frozen reference-side separating inequalities for crossing triangle pairs.
This is an experimental untangling search, not friction, tissue contact physics
or continuous collision detection. Never changes mesh connectivity or pins.
Each pair uses a separating axis from a screened reference pose, transported
by a proper rigid fit of the six points. All 9 cross-vertex inequalities
are imposed; satisfying them implies separation along that frozen axis.
"""
import numpy as np
from scipy.sparse import coo_matrix,diags

def reference_axis(A,B):
    ea=np.roll(A,-1,axis=0)-A;eb=np.roll(B,-1,axis=0)-B
    na=np.cross(ea[0],ea[1]);nb=np.cross(eb[0],eb[1])
    axes=np.vstack([na,nb,np.cross(ea[:,None,:],eb[None,:,:]).reshape(-1,3),
                    np.cross(na,ea),np.cross(nb,eb)])
    norms=np.linalg.norm(axes,axis=1);axes=axes[norms>1e-14]/norms[norms>1e-14,None]
    axes=np.vstack([axes,-axes])
    gaps=(B@axes.T).min(0)-(A@axes.T).max(0)
    i=int(gaps.argmax())
    if gaps[i]<=1e-10:raise ValueError('Reference pair has no verified positive separating axis')
    return axes[i],float(gaps[i])

def transported_axis(A,B,C,D):
    n,gap=reference_axis(A,B)
    p=np.vstack([A,B]);q=np.vstack([C,D])
    U,_,Vt=np.linalg.svd((p-p.mean(0)).T@(q-q.mean(0)))
    R=U@np.diag([1,1,np.linalg.det(U@Vt)])@Vt
    return n@R,gap

class SelfSeparation:
    def __init__(self,solid,faces,reference,current,pairs,rho):
        assert rho>0 and np.isfinite(rho)
        self.shape=current.shape;self.rho=float(rho)
        rr=[];cc=[];vv=[];margin=[];self.details=[]
        for i,j in sorted(pairs):
            a,b=faces[i],faces[j]
            assert not np.intersect1d(a,b).size
            n,gap=transported_axis(reference[a],reference[b],current[a],current[b])
            m=min(.02,.2*gap)
            self.details.append({'pair':[int(i),int(j)],'referenceGapMm':gap,'marginMm':m,'normal':n.tolist()})
            for u in a:
                for v in b:
                    k=len(margin);margin.append(m)
                    rr.extend([k]*6);cc.extend(np.r_[u*3+np.arange(3),v*3+np.arange(3)]);vv.extend(np.r_[-n,n])
        self.A=coo_matrix((vv,(rr,cc)),shape=(len(margin),current.size)).tocsr()
        self.margin=np.asarray(margin)
        self.weights=np.full(len(margin),1/max(len(margin),1))
        free=(solid.free[:,None]*3+np.arange(3)).ravel()
        fixed_rows=np.diff(self.A[:,free].tocsr().indptr)==0
        self.fixedConflictCount=int((fixed_rows&(self.A@current.ravel()<self.margin-1e-9)).sum())
    def violation(self,x):
        return np.maximum(self.margin-self.A@x.ravel(),0.)
    def evaluate(self,x,gradient=True):
        v=self.violation(x)
        value=.5*self.rho*np.dot(self.weights,v*v)
        g=-(self.A.T@(self.rho*self.weights*v)).reshape(self.shape) if gradient else None
        return float(value),g
    def hessian(self,x):
        return self.A.T@diags(self.rho*self.weights*(self.violation(x)>0))@self.A
