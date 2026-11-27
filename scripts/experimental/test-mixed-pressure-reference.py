import numpy as np
from volumetric_reference import Solid
from mixed_pressure_reference import MixedPressure

rest=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,-1]])
solid=Solid(rest,np.array([[0,1,2,3],[0,2,1,4]]),np.array([0,3,4]))
x=rest.copy();x[3,2]=1.2;x[4,2]=-.85
zero=MixedPressure(solid,jump=0)
E,g,J=zero.evaluate(x,30)
E0,g0,J0=solid.evaluate(x,90)
assert abs(E-E0)<1e-12 and np.max(np.abs(g-g0))<1e-11
mix=MixedPressure(solid,jump=1)
assert np.allclose(mix.pressure(np.ones(2)*.2)[0],12)
assert np.allclose(mix.pressure(np.array([.2,-.2]))[0],[4,-4])
for k in range(x.size):
    plus,minus=x.copy(),x.copy();plus.flat[k]+=1e-6;minus.flat[k]-=1e-6
    finite=(mix.evaluate(plus,30,False)[0]-mix.evaluate(minus,30,False)[0])/2e-6
    assert abs(finite-mix.evaluate(x,30)[1].flat[k])<1e-7
sol,report=mix.minimize(x,iterations=80)
assert report['converged'],report
assert np.array_equal(sol[solid.fixed],x[solid.fixed])
assert mix.audit(sol)['JRange'][0]>0 and mix.audit(sol)['pressureMeanBalanceError']<1e-12
failed,report=mix.minimize(x,iterations=0,tolerance_mm=1e-14)
assert not report['converged'] and report['rolledBack'] and np.array_equal(failed,x)
bad=Solid(rest,solid.tets,solid.fixed,np.array([.7,.3]))
try:MixedPressure(bad)
except ValueError:pass
else:raise AssertionError('Nonuniform integration accepted')
print('PASS mixed pressure: zero-jump exact bulk, constant/opposite pressure oracle, exact reduced derivative, pins, positive J, rollback, cut-cell refusal')
