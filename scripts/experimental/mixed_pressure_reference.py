"""Experimental P1 displacement / P0 pressure with face-jump stabilization.

Inspired by Becker, Burman and Hansbo (2009), https://research.chalmers.se/en/publication/97621.
This is a finite-deformation log-J pressure experiment, NOT their linear
elasticity formulation or a proof of mesh-independent inf-sup stability.
A retained local bulk term prevents a purely averaged volume objective.
Only fully physical tetrahedra are supported; cut cells need new integration.
"""
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import splu
from interior_face_refinement import face_owners
from local_volume_constraints import log_volume_jacobian
from projected_solid_metric import factor


class MixedPressure:
    def __init__(self, solid, bulk_pressure=60., jump=1.):
        self.solid=solid
        self.bulk_pressure=float(bulk_pressure)
        self.jump=float(jump)
        if not np.isfinite([self.bulk_pressure,self.jump]).all() or self.bulk_pressure<=0 or self.jump<0:
            raise ValueError('Invalid pressure parameters')
        if not np.allclose(solid.weights,solid.volumes/solid.volumes.sum(),rtol=1e-11,atol=1e-15):
            raise ValueError('Mixed reference requires full physical tetrahedra')
        self.cells=np.arange(len(solid.tets))
        pairs=np.asarray([v for v in face_owners(solid.tets).values() if len(v)==2],dtype=int).reshape(-1,2)
        self.pairs=pairs
        a,b=pairs.T
        w=solid.weights
        # Harmonic incident volume = area times harmonic normal height / 3,
        # normalized by total rest volume. Symmetric, positive, no muscle coefficients.
        self.face_weights=2*w[a]*w[b]/(w[a]+w[b])
        k=self.face_weights*self.jump/self.bulk_pressure
        self.S=coo_matrix((np.r_[k,k,-k,-k],(np.r_[a,b,a,b],np.r_[a,b,b,a])),shape=(len(w),len(w))).tocsc()
        self.M=(diags(w/self.bulk_pressure)+self.S).tocsc()
        self.linear=splu(self.M,permc_spec='MMD_AT_PLUS_A')
        self.matrix_norm=float(np.asarray(abs(self.M).sum(axis=1)).max())

    def pressure(self, logj):
        rhs=self.solid.weights*logj
        p=self.linear.solve(rhs)
        residual=self.M@p-rhs
        scale=self.matrix_norm*np.max(np.abs(p))+np.max(np.abs(rhs))
        error=float(np.max(np.abs(residual))/max(scale,1e-300))
        if not np.isfinite(p).all() or not np.isfinite(error) or error>1e-10:
            raise RuntimeError('Mixed pressure solve failed backward error')
        return p,error

    def evaluate(self,x,bulk=30.,gradient=True):
        # E_iso + local bulk + max_p(p.W.logJ - 1/2 p.M.p).
        energy,g,J=self.solid.evaluate(x,bulk,gradient)
        if not np.isfinite(energy):return energy,g,J
        c=np.log(J)
        p,_=self.pressure(c)
        energy+=.5*np.dot(self.solid.weights*c,p)
        if gradient:
            _,A=log_volume_jacobian(self.solid,x,self.cells)
            g+=(A.T@(self.solid.weights*p)).reshape(x.shape)
        return float(energy),g,J

    def audit(self,x,bulk=30.):
        J=np.linalg.det(self.solid.edges(x)@self.solid.inverse)
        if np.any(J<=1e-9) or not np.isfinite(J).all():raise ValueError('Invalid audit J')
        c=np.log(J);p,error=self.pressure(c)
        a,b=self.pairs.T
        def jump_rms(q):
            return float(np.sqrt(np.dot(self.face_weights,(q[a]-q[b])**2)/self.face_weights.sum())) if len(a) else 0.
        return {'JRange':[float(J.min()),float(J.max())],
                'maxAbsLogJ':float(np.max(np.abs(c))),
                'globalVolumeRatio':float(np.dot(self.solid.weights,J)),
                'pressureRange':[float(p.min()),float(p.max())],
                'pressureJumpRms':jump_rms(p),
                'unfilteredPressureJumpRms':jump_rms(self.bulk_pressure*c),
                'pressureBackwardError':error,
                'pressureMeanBalanceError':float(abs(np.dot(self.solid.weights,p)/self.bulk_pressure-np.dot(self.solid.weights,c))),
                'pressureDofs':len(p),'interiorFaces':len(a),'infSupCertified':False,
                'allCellsAudited':True,'localBulk':bulk,'pressureBulk':self.bulk_pressure,'jump':self.jump}

    def minimize(self,initial,bulk=30.,iterations=40,tolerance_mm=.001):
        """Exact reduced objective, local bulk tangent as search metric only.
        Pressure is recomputed on every trial; no frozen pressure acceptance.
        Nonconvergence returns original input and retains trial diagnostics.
        """
        x=initial.copy();s=self.solid;rows=[];reason='iteration-limit'
        for i in range(iterations):
            E,g,J=self.evaluate(x,bulk)
            if not np.isfinite(E):raise ValueError('Nonpositive initial state')
            residual=float(np.linalg.norm(s.factor.solve(g[s.free]),axis=1).max(initial=0))
            if residual<tolerance_mm:reason='residual-tolerance';break
            success=False
            for project in (False,True):
                try:
                    lin,_=factor(s,x,bulk+self.bulk_pressure,project)
                    d=lin.solve(-g[s.free].ravel()).reshape(-1,3)
                except RuntimeError:continue
                slope=float(np.sum(g[s.free]*d))
                if not np.isfinite(d).all() or slope>=0:continue
                gain=1.
                for _ in range(40):
                    trial=x.copy();trial[s.free]+=gain*d
                    val,_,_=self.evaluate(trial,bulk,False)
                    if np.isfinite(val) and val<=E+.0001*gain*slope:
                        success=True;break
                    gain*=.5
                if success:break
            if not success:reason='line-search-stalled';break
            x=trial
            rows.append({'iteration':i+1,'energy':val,'gain':gain,'previousResidualMm':residual})
        E,g,J=self.evaluate(x,bulk)
        residual=float(np.linalg.norm(s.factor.solve(g[s.free]),axis=1).max(initial=0))
        if residual<tolerance_mm:reason='residual-tolerance'
        ok=reason=='residual-tolerance'
        assert np.array_equal(x[s.fixed],initial[s.fixed])
        return (x if ok else initial.copy()),{'converged':ok,'reason':reason,'residualMm':residual,
            'iterations':len(rows),'trialAudit':self.audit(x,bulk),'rows':rows,
            'rolledBack':not ok,'accepted':False,'boneContactEvaluated':False}
