"""Offline local source-triangle shoulder optimization, not a runtime solver.

A nearest-point linearization supplies a search direction, not a collision
certificate. Exact all-source vertex and both edge directions must screen
each resulting trial independently. No joint translations or altered bone size.
"""
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation as R
from scipy.optimize import minimize
import trimesh


class LinkedShoulder:
    def __init__(self,data):
        f=data['frames'];self.L={k:np.asarray(v)/data['calibration']['scale'] for k,v in data['calibration']['landmarks'].items()}
        self.T,self.C,self.S=(R.from_quat(f[k]) for k in ('thorax','clavicle','scapula'))
        self.sc=np.asarray(f['SC0']);self.ac=np.asarray(f['AC0'])
    def forward(self,a):
        a=np.deg2rad(a)
        C=self.T*R.from_euler('YXZ',self.sc+a[:3]*[-1,-1,1])
        S=C*R.from_euler('YXZ',self.ac+a[3:]*[1,-1,1])
        qc=C*self.C.inv();qs=S*self.S.inv()
        AC=self.L['SC']+qc.apply(self.L['AC']-self.L['SC'])
        GH=AC+qs.apply(self.L['GH']-self.L['AC'])
        return qc,qs,AC,GH
    def pose_points(self,a,p):
        _,q,AC,_=self.forward(a)
        return q.apply(p-self.L['AC'])+AC
    def verify(self,row,scale):
        qc,qs,AC,GH=self.forward(row['parameters'])
        for name,q,p in [('clavicle',qc,self.L['SC']),('scapula',qs,AC),('humerus',None,GH)]:
            b=row['bones'][name]
            assert np.linalg.norm(np.asarray(b['p'])/scale-p)<1e-7
            if q is not None:
                assert (q.inv()*R.from_quat(b['q'])).magnitude()<1e-9
    def candidate(self,row,delta,scale):
        a=np.asarray(row['parameters'])+delta
        qc,qs,AC,GH=self.forward(a)
        bones={**row['bones']}
        for name,q,p in [('clavicle',qc,self.L['SC']),('scapula',qs,AC),('humerus',None,GH)]:
            b=dict(bones[name]);b['p']=(p*scale).tolist()
            if q is not None:b['q']=q.as_quat().tolist()
            bones[name]=b
        # All mesh origins in exported linked hierarchy are joint pivots.
        assert np.linalg.norm(np.asarray(bones['scapula']['rest'])/scale-self.L['AC'])<1e-7
        assert np.linalg.norm(np.asarray(bones['humerus']['rest'])/scale-self.L['GH'])<1e-7
        return {**row,'parameters':a.tolist(),'bones':bones,
                'priorDeviationL2Deg':float(np.linalg.norm(a-np.asarray(row['prior'])))}


def nearest(mesh,p):
    results=[trimesh.proximity.closest_point(mesh,p[i:i+128]) for i in range(0,len(p),128)]
    return np.concatenate([r[0] for r in results]),np.concatenate([r[1] for r in results])


def direction(rig,row,scap,thorax,costal,rest_gap,trust=1.):
    a=np.asarray(row['parameters']);prior=np.asarray(row['prior'])
    # Use actual closest triangles for residuals. Surface-vertex KD search
    # selects near-contact candidates only; it is not the final collision test.
    rest=np.asarray(scap.vertices);v=rig.pose_points(a,rest)
    rv=np.asarray(thorax.vertices);tree=cKDTree(rv)
    d,_=tree.query(v)
    ids=np.argsort(d,kind='stable')[:min(512,len(v))]
    _,reverse=cKDTree(v).query(rv)
    # Include one representative from each closest scapular vertex to ribs;
    # limit to 512 most proximal additional vertices deterministically.
    reverse=np.unique(reverse)
    reverse=reverse[np.argsort(d[reverse],kind='stable')[:512]]
    ids=np.unique(np.r_[ids,reverse])
    p=np.vstack([costal,rest[ids]])
    x=rig.pose_points(a,p);q,gap=nearest(thorax,x)
    if np.any(gap<1e-8):return None,{'reason':'ambiguous-zero-nearest-vector'}
    n=(x-q)/gap[:,None]
    # Derivative of true distance away from nearest-feature switching.
    jac=np.empty((len(p),6))
    for j in range(6):
        da=np.zeros(6);da[j]=1e-4
        dx=(rig.pose_points(a+da,p)-rig.pose_points(a-da,p))/(2e-4)
        jac[:,j]=np.einsum('ij,ij->i',n,dx)
    count=len(costal);e=gap[:count]-rest_gap;A=jac[:count];B=jac[count:];g=gap[count:]
    # No greater unweighted departure from the original literature prior.
    radius=np.linalg.norm(a-prior)
    objective=lambda u:float(np.mean((e+A@u)**2))
    gradient=lambda u:2*A.T@(e+A@u)/count
    constraints=[{'type':'ineq','fun':lambda u:g+B@u,'jac':lambda u:B},
                 {'type':'ineq','fun':lambda u:radius**2-np.dot(a+u-prior,a+u-prior),
                  'jac':lambda u:-2*(a+u-prior)}]
    z=minimize(objective,np.zeros(6),jac=gradient,bounds=[(-trust,trust)]*6,
               constraints=constraints,method='SLSQP',options={'maxiter':80,'ftol':1e-11})
    info={'trustDeg':trust,'localSampleCount':len(ids),'costalSampleCount':count,
          'linearizedObjectiveBefore':objective(np.zeros(6)),
          'linearizedObjectiveAfter':objective(z.x),'solverSuccess':bool(z.success),
          'solverMessage':str(z.message),'minimumLinearizedGapMm':float(np.min(g+B@z.x)),
          'priorDeviationAfterDeg':float(np.linalg.norm(a+z.x-prior)),
          'certificate':False}
    if not z.success or np.min(g+B@z.x)<-1e-7 or np.linalg.norm(a+z.x-prior)>radius+1e-7 or objective(z.x)>=objective(np.zeros(6))-1e-10:
        return None,info
    return z.x,info
