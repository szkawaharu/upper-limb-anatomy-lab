"""Check interior interpolation independently with a known affine field."""
import importlib.util
from pathlib import Path
import numpy as np
spec=importlib.util.spec_from_file_location('transfer',Path(__file__).with_name('transfer-volumetric-reference.py'));module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
p=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1.]],float);cells=np.array([[0,1,2,3]]);new=np.vstack([p,[.2,.3,.1],[2,2,2]]);new_cells=np.array([[0,1,2,4],[0,1,4,3],[0,4,2,3],[4,1,2,3]])
A=np.array([[1.1,.2,.1],[.1,.9,.2],[0,.1,1.05]]);at=np.array([.3,-.4,.5]);current=p@A.T+at
out,report=module.transfer(p,cells,current,new,4,new_cells)
assert np.array_equal(out[:4],current)and np.linalg.norm(out[4]-(new[4]@A.T+at))<1e-12
assert np.array_equal(out[5],new[5])and report['transferredInteriorNodes']==1
bad=new.copy();bad[4]=[2,2,2]
try:module.transfer(p,cells,current,bad,4,new_cells)
except ValueError:pass
else:raise AssertionError('An outside node was snapped to a nearest element')
bad=new.copy();bad[0,0]+=1e-8
try:module.transfer(p,cells,current,bad,4,new_cells)
except ValueError:pass
else:raise AssertionError('Changed anatomical source vertex was accepted')
print('PASS transfer: independent affine reproduction, exact source boundary, inactive generated node retained, outside points and changed source rejected')
