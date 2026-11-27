import numpy as np
from continuation import advance
attempts=[]
def solve(x,previous,target):
 assert abs(x[0]-previous)<1e-14,'Unconverged state was reused'
 attempts.append((previous,target))
 if target-previous>.15:
  x[:]=999;return x,{'converged':False,'reason':'synthetic-iteration-limit','residualMaxMm':2.}
 x[0]=target;return x,{'converged':True,'residualMaxMm':.0005}
x,r=advance([0.],solve,.4);assert r['completed']and x[0]==1. and len(r['rejected'])>0
assert all(s['residualMaxMm']<.001 for s in r['accepted'])
def impossible(x,previous,target):return x+100,{'converged':False,'reason':'no-solution'}
x,r=advance([0.],impossible,.5,max_rejections=3);assert not r['completed']and r['progress']==0 and x[0]==0 and len(r['rejected'])==3
print('PASS numerical continuation: rollback unconverged trial, preserve accepted boundary, exact endpoint, fail without false progress')
