"""Affine-invariant minimum-bending extension for unsupported spline controls.
It supplies no anatomical mass. Condensation and numerical stabilization are
explicitly reported; material integration must be independently validated.
"""
import itertools
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

class BendingExtension:
 def __init__(self,cage,active):
  self.active=np.unique(active);self.inactive=np.setdiff1d(np.arange(len(cage.nodes)),self.active);rows=[];cols=[];values=[];row=0
  shape=cage.shape
  for axis in range(3):
   ranges=[range(int(n))for n in shape];ranges[axis]=range(1,int(shape[axis])-1)
   for i in itertools.product(*ranges):
    for offset,coef in [(-1,1),(0,-2),(1,1)]:
     p=list(i);p[axis]+=offset;rows.append(row);cols.append(np.ravel_multi_index(p,shape));values.append(coef/cage.h)
    row+=1
  for a,b in itertools.combinations(range(3),2):
   ranges=[range(int(n))for n in shape];ranges[a]=range(int(shape[a])-1);ranges[b]=range(int(shape[b])-1)
   for i in itertools.product(*ranges):
    for da,db,coef in [(0,0,1),(1,0,-1),(0,1,-1),(1,1,1)]:
     p=list(i);p[a]+=da;p[b]+=db;rows.append(row);cols.append(np.ravel_multi_index(p,shape));values.append(np.sqrt(2)*coef/cage.h)
    row+=1
  B=coo_matrix((values,(rows,cols)),shape=(row,len(cage.nodes))).tocsr();self.full=(B.T@B/row).tocsc();A=self.active;I=self.inactive
  self.cross=self.full[I][:,A];self.factor=splu(self.full[I][:,I])if len(I)else None
  self.extension=-self.factor.solve(self.cross.toarray())if len(I)else np.zeros((0,len(A)))
  self.condensed=self.full[A][:,A].toarray()+self.full[A][:,I]@self.extension
  self.condensed=(self.condensed+self.condensed.T)/2
  self.index=np.full(len(cage.nodes),-1,int);self.index[A]=np.arange(len(A))

 def extend(self,x,rest):
  out=x.copy();out[self.inactive]=rest[self.inactive]+self.extension@(x[self.active]-rest[self.active]);return out
