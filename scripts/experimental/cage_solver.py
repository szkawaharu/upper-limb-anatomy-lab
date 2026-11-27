"""Direct sparse SPD quasi-Newton reference for a small embedding cage.
Exact normalized material energy; positive-J Armijo acceptance. Offline only.
"""
import numpy as np
from projected_solid_metric import factor

def minimize(solid,initial,bulk=30.,iterations=80,tolerance_mm=.001,observer=None):
 x=initial.copy();rows=[];reason='iteration-limit'
 for k in range(iterations):
  energy,g,J=solid.evaluate(x,bulk)
  if not np.isfinite(energy):raise ValueError('Non-positive input cage')
  residual=float(np.linalg.norm(solid.factor.solve(g[solid.free]),axis=1).max(initial=0))
  if residual<tolerance_mm:reason='residual-tolerance';break
  success=False;trial=x.copy()
  # A constrained global tangent can be positive even when individual cell
  # tangents are indefinite. Try the exact Newton direction first; retain
  # projected-local fallback and the identical energy/J/Armijo safeguards.
  for project in [False,True]:
   try:
    linear_factor,metric_info=factor(solid,x,bulk,project);d=linear_factor.solve(-g[solid.free].ravel()).reshape(-1,3)
   except RuntimeError:continue
   slope=float(np.sum(g[solid.free]*d));gain=1.
   if not np.isfinite(d).all()or not slope<0:continue
   for _ in range(40):
    trial[solid.free]=x[solid.free]+gain*d;value,_,nextJ=solid.evaluate(trial,bulk,gradient=False)
    if np.isfinite(value)and value<=energy+.0001*gain*slope:success=True;break
    gain*=.5
   if success:break
  if not success:reason='line-search-stalled';break
  x=trial;rows.append({'iteration':k+1,'energy':value,'gain':gain,'minJ':float(nextJ.min()),'maxJ':float(nextJ.max()),'previousResidualMm':residual,'projectedMetric':project})
  if observer:observer(x,rows[-1])
 energy,g,J=solid.evaluate(x,bulk);residual=float(np.linalg.norm(solid.factor.solve(g[solid.free]),axis=1).max(initial=0))
 if residual<tolerance_mm:reason='residual-tolerance'
 return x,{'converged':reason=='residual-tolerance','reason':reason,'iterations':len(rows),'residualMaxMm':residual,'energy':energy,'cageJQuantiles':np.quantile(J,[0,.01,.5,.99,1]).tolist(),'cageNonPositive':int((J<=0).sum()),'quadratureVolumeRatio':float(np.dot(solid.weights,J)),'rows':rows,'metric':locals().get('metric_info')}
