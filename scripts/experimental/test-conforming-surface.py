import numpy as np
from tetrahedral_cage import GridCage
from conforming_surface import conform
p=np.array([[.17,.13,.23],[3.37,.22,.41],[.32,2.91,.49],[.22,.34,3.15]])
f=np.array([[0,2,1],[0,1,3],[1,2,3],[2,0,3]])
c=GridCage(p,.8);r=conform(c,p,f);v=r['points'];tri=r['faces'];assert np.array_equal(v[:len(p)],p)
assert r['sourceAreaMaxRelativeError']<1e-10
edges=np.sort(np.concatenate([tri[:,[0,1]],tri[:,[1,2]],tri[:,[2,0]]]),axis=1);_,count=np.unique(edges,axis=0,return_counts=True);assert np.all(count==2)
ids,bary,_=c.embed(v);nodes=c.nodes+np.random.default_rng(14).normal(0,.015,c.nodes.shape)
E=lambda x:np.stack([x[c.tets[:,i]]-x[c.tets[:,0]]for i in [1,2,3]],axis=2)
F=E(nodes)@np.linalg.inv(E(c.nodes));assert np.linalg.det(F).min()>.5
now=c.positions(nodes,ids,bary);normal=lambda x:np.cross(x[tri[:,1]]-x[tri[:,0]],x[tri[:,2]]-x[tri[:,0]])
field=F[r['element']];expected=np.linalg.det(field)[:,None]*np.einsum('tij,tj->ti',np.linalg.inv(field).transpose(0,2,1),normal(v))
assert np.max(np.linalg.norm(normal(now)-expected,axis=1))<1e-12
assert np.all(np.einsum('ti,ti->t',normal(now),expected)>0)
print('PASS conforming source: all original vertices retained, per-source triangle area preserved, closed oriented subdivision, exact affine-cell surface normals under nonuniform positive-J warp')
