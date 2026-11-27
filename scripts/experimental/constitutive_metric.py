"""Exact derivative of the normalized isochoric + log-volume stress.
The optimization metric projects negative local eigenvalues to zero; objective,
stress, source constraints and positive-J acceptance are never changed.
"""
import numpy as np

def tangent(F,bulk=30.,project=True):
 J=np.linalg.det(F)
 if J.min()<=0:raise ValueError('Non-positive constitutive input')
 C=np.linalg.inv(F).transpose(0,2,1);Ic=np.einsum('qij,qij->q',F,F);a=J**(-2/3);f=F.reshape(-1,9);c=C.reshape(-1,9)
 H=a[:,None,None]*np.eye(9)-(2/3)*a[:,None,None]*(f[:,:,None]*c[:,None,:]+c[:,:,None]*f[:,None,:])
 H+=(bulk+2*a*Ic/9)[:,None,None]*c[:,:,None]*c[:,None,:]
 T=np.einsum('qil,qkj->qijkl',C,C).reshape(-1,9,9)
 H-=(bulk*np.log(J)-a*Ic/3)[:,None,None]*T;H=(H+H.transpose(0,2,1))/2
 if not project:return H
 eigen,U=np.linalg.eigh(H);root=U*np.sqrt(np.maximum(eigen,0))[:,None,:]
 return root,{'negativeEigenvalues':int((eigen<0).sum()),'minUnprojectedEigenvalue':float(eigen.min()),'projection':'negative local material eigenvalues to zero; exact energy and stress retained'}
