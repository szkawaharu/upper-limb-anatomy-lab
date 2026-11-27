"""Adaptive numerical loading with transactional accepted states.
An unconverged nonlinear trial never becomes the next initial state. Numerical
substeps do not change the endpoint, objective, boundary data or tolerance.
"""
import numpy as np

def advance(initial,solve,base_step,on_accept=None,min_step=1e-7,max_rejections=16):
 x=np.array(initial,copy=True);progress=0.;step=base_step;accepted=[];rejected=[];consecutive=0
 while progress<1-1e-12:
  target=min(1.,progress+step);trial,report=solve(x.copy(),progress,target)
  if report.get('converged')and np.isfinite(trial).all():
   x=np.array(trial,copy=True);progress=target;accepted.append({'progress':progress,**report});consecutive=0
   if on_accept:on_accept(x,progress,accepted,rejected)
   step=min(base_step,step*1.25)
  else:
   rejected.append({'fromProgress':progress,'attemptedProgress':target,**report});consecutive+=1;step*=.5
   if step<min_step or consecutive>=max_rejections:break
 return x,{'progress':progress,'completed':progress>=1-1e-12,'accepted':accepted,'rejected':rejected}
