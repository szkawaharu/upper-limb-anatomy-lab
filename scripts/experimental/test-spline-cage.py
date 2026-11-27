import numpy as np
from spline_cage import SplineCage,SplineSolid
rng=np.random.default_rng(21);points=rng.uniform(-2,2,(40,3));c=SplineCage(points,.8);embedding=c.embed(points);ids,w,g=embedding
assert np.max(np.abs(w.sum(1)-1))<1e-13 and np.max(np.abs(g.sum(1)))<1e-13
assert np.max(np.abs(c.positions(c.nodes,embedding)-points))<1e-12
A=np.array([[1.2,.1,.2],[0,.8,.1],[.2,.1,1.1]]);b=np.array([2,-3,4]);x=c.nodes@A.T+b
assert np.max(np.abs(c.positions(x,embedding)-(points@A.T+b)))<1e-12
assert np.max(np.abs(c.gradients(x,embedding)-A))<1e-12
# Independently finite-difference the spatial derivative across support changes.
x=c.nodes+rng.normal(0,.01,c.nodes.shape)
for k in range(3):
 delta=np.eye(3)[k]*1e-6;fd=(c.positions(x,c.embed(points+delta))-c.positions(x,c.embed(points-delta)))/2e-6
 assert np.max(np.abs(fd-c.gradients(x,embedding)[:,:,k]))<1e-8
for location in [.5,1.5,2.5]:
 p=np.array([[location*.8,0,0]]);left=c.gradients(x,c.embed(p-[1e-8,0,0]));right=c.gradients(x,c.embed(p+[1e-8,0,0]));assert np.max(np.abs(left-right))<1e-7
# All occupied quadrature nodes fixed except a well-conditioned interior block.
classes=np.zeros(len(c.nodes),np.int8);classes[((c.nodes>-.4)&(c.nodes<.4)).all(1)]=-1
s=SplineSolid(c,points,np.ones(len(points)),classes);E,grad,J,H=s.evaluate(x,30,True);v=rng.normal(size=x.shape);eps=1e-6
fd=(s.evaluate(x+eps*v,30)[0]-s.evaluate(x-eps*v,30)[0])/2/eps
assert abs(fd-np.sum(grad*v))<1e-7
assert np.max(np.abs(H-H.T),initial=0)<1e-12
print('PASS C1 spline: source/affine reproduction, exact spatial gradient, support-change continuity, independent energy gradient, symmetric strain metric')
