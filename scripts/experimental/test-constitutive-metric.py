import numpy as np
from constitutive_metric import tangent
rng=np.random.default_rng(9);F=np.eye(3)[None]+rng.normal(0,.1,(7,3,3));H=tangent(F,30,False)
def stress(F):
 J=np.linalg.det(F);C=np.linalg.inv(F).transpose(0,2,1);Ic=np.einsum('qij,qij->q',F,F)
 return J[:,None,None]**(-2/3)*(F-Ic[:,None,None]*C/3)+30*np.log(J)[:,None,None]*C
for k in range(9):
 e=np.zeros((3,3));e.flat[k]=1e-6;fd=((stress(F+e)-stress(F-e))/2e-6).reshape(-1,9)
 assert np.max(np.abs(fd-H[:,:,k]))<1e-7
assert np.max(np.abs(H-H.transpose(0,2,1)))<1e-12
root,report=tangent(F,30,True);P=root@root.transpose(0,2,1);assert np.linalg.eigvalsh(P).min()>-1e-12
# Independent rest tangent, including exact infinitesimal rotational null modes.
R=tangent(np.eye(3)[None],30,False)[0];W=np.array([[0,-.3,.2],[.3,0,-.1],[-.2,.1,0]])
assert np.linalg.norm(R@W.ravel())<1e-12
print('PASS material tangent: independent stress derivative, symmetry, positive projection, exact rest rotational null modes')
