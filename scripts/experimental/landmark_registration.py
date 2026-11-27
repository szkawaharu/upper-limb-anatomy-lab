"""Proper, non-reflecting landmark registration for cross-atlas diagnostics.
Uniform scales do not morph anatomy. Residuals and uncertain anchors must remain
visible; registered point paths do not define broad enthesis footprints.
"""
import numpy as np

def similarity(source,target):
 a,b=np.asarray(source,float),np.asarray(target,float)
 if a.shape!=b.shape or a.ndim!=2 or a.shape[1]!=3 or len(a)<3:raise ValueError('Three matching 3-D anchors required')
 if not np.isfinite(a).all()or not np.isfinite(b).all():raise ValueError('Nonfinite landmarks')
 ac,bc=a.mean(0),b.mean(0);A,B=a-ac,b-bc
 if np.linalg.matrix_rank(A)<2 or np.linalg.matrix_rank(B)<2:raise ValueError('Collinear landmarks')
 U,s,Vt=np.linalg.svd(A.T@B);sign=np.linalg.det(Vt.T@U.T);D=np.diag([1,1,sign]);R=Vt.T@D@U.T
 scale=float(np.sum(s*np.diag(D))/np.sum(A*A))
 if not scale>0:raise ValueError('Nonpositive registration scale')
 t=bc-scale*(R@ac);fit=scale*a@R.T+t
 return {'rotation':R,'scale':scale,'translation':t,'errors':np.linalg.norm(fit-b,axis=1)}

def bone_similarity(source_prox,source_dist,source_lateral,target_prox,target_dist,target_lateral):
 def frame(prox,dist,lateral):
  y=np.asarray(prox)-dist;length=np.linalg.norm(y)
  if length<1e-12:raise ValueError('Zero bone length')
  y=y/length;x=np.cross(y,lateral);n=np.linalg.norm(x)
  if n<1e-12:raise ValueError('Degenerate transverse bone direction')
  x=x/n;z=np.cross(x,y);return np.column_stack([x,y,z]),length
 A,la=frame(source_prox,source_dist,source_lateral);B,lb=frame(target_prox,target_dist,target_lateral)
 R=B@A.T;scale=lb/la;t=np.asarray(target_prox)-scale*R@source_prox
 return {'rotation':R,'scale':scale,'translation':t}

def apply(registration,points):
 return registration['scale']*np.asarray(points)@registration['rotation'].T+registration['translation']
