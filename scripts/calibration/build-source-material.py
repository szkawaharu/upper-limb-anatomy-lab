"""Precompute static rest centres, smooth source half-chord vector fields and
closed coarse volume quadrature meshes. No animation-dependent rest fitting.
"""
import json,sys
from pathlib import Path
import numpy as np,trimesh
from scipy.sparse import coo_matrix,diags
src=Path(sys.argv[1]if len(sys.argv)>1 else'test-results/v076-source-clean/bindings.json');out=Path(sys.argv[2]if len(sys.argv)>2 else'test-results/v076-source-material');out.mkdir(parents=True,exist_ok=True)
records=[]
for B in json.loads(src.read_text()):
 id=B['spec']['id'];p=np.array(B['rest']).reshape(-1,3);faces=np.array(B['faces']).reshape(-1,3);n=len(p)
 N=np.zeros_like(p);fn=np.cross(p[faces[:,1]]-p[faces[:,0]],p[faces[:,2]]-p[faces[:,0]])
 for k in range(3):np.add.at(N,faces[:,k],fn)
 N/=np.maximum(np.linalg.norm(N,axis=1,keepdims=True),1e-20)
 h=np.fromfile('test-results/v076-source-clean/thickness/'+id+'.f32',dtype='<f4');assert len(h)==n
 d0=N*h[:,None];D=d0.copy()
 edges=np.unique(np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0);a,b=edges.T;weights=1/np.maximum(np.linalg.norm(p[a]-p[b],axis=1),1e-6)
 W=coo_matrix((np.r_[weights,weights],(np.r_[a,b],np.r_[b,a])),shape=(n,n)).tocsr();den=np.asarray(W.sum(1)).ravel();W=diags(1/np.maximum(den,1e-20))@W
 for _ in range(80):D=.99*(W@D)+.01*d0
 controls=np.array(B['field']['controls']);t=np.array(B['material']);l=np.clip(B['lane'],0,len(controls)-1);i=np.floor(l).astype(int);u=l-i;v=1-u
 laneW=np.array([v**3/6,(3*u**3-6*u*u+4)/6,(-3*u**3+3*u*u+3*u+1)/6,u**3/6]);bezier=np.array([(1-t)**3,3*t*(1-t)**2,3*t*t*(1-t),t**3]);C=np.zeros_like(p)
 for k in range(4):C+=np.einsum('nij,ni->nj',controls[np.clip(i+k-1,0,len(controls)-1)],bezier.T)*laneW[k,:,None]
 data=np.concatenate([C.ravel(),D.ravel(),N.ravel()]).astype('<f4');(out/(id+'.f32')).write_bytes(data.tobytes())
 mesh=trimesh.Trimesh(vertices=p,faces=faces,process=False)
 quadrature=None
 if mesh.is_watertight and mesh.is_winding_consistent:
  w=np.array(B['weights']);variable=B['spec']['id']in['v58Levator','v59TrapUpper'];pureA=(w[faces]==0).all(1)if not variable else np.zeros(len(faces),bool);pureB=(w[faces]==1).all(1);mixed=~(pureA|pureB)
  area=np.linalg.norm(fn,axis=1);centroids=p[faces].mean(1);sourceCenter=p.mean(0);target=500;cell=np.sqrt(area[mixed].sum()/2/target)
  axis=np.argmax(np.abs(fn),axis=1);orient=axis*2+(fn[np.arange(len(fn)),axis]>0);spatial=np.floor((centroids-p.min(0))/max(cell,1e-6)).astype(int);faceW=w[faces].mean(1);bins=np.floor(faceW*8).astype(int)
  groups={}
  for i in np.flatnonzero(mixed & (area>1e-15)):
   key=(*spatial[i],orient[i],bins[i]);groups.setdefault(key,[]).append(i)
  chosen=[];multipliers=[]
  for ids in groups.values():
   ids=np.array(ids);mean=np.average(centroids[ids],axis=0,weights=area[ids]);j=ids[np.argmin(np.linalg.norm(centroids[ids]-mean,axis=1))];chosen.append(j);multipliers.append(float(area[ids].sum()/area[j]))
  chosen=np.array(chosen);qfaces=faces[chosen];used=np.unique(qfaces);remap=np.full(n,-1,int);remap[used]=np.arange(len(used));qfaces=remap[qfaces];R=p[faces[chosen]]-sourceCenter;control=np.einsum('ij,ij->i',R[:,0],np.cross(R[:,1],R[:,2]))
  quadrature={'sourceVertices':used.tolist(),'faces':qfaces.ravel().tolist(),'multipliers':multipliers,'controlDet':control.tolist(),'referenceVolume':float(mesh.volume),'sourceCenter':sourceCenter.tolist(),'pureBVector':fn[pureB].sum(0).tolist(),'ignoredRigidOriginFaces':int(pureA.sum()),'exactRigidInsertionFaces':int(pureB.sum()),'method':'stratified source-surface integral of volume increments relative to a rigid control variate','cellWorld':float(cell),'nominalBudget':target}
 (out/(id+'-quadrature.json')).write_text(json.dumps(quadrature,separators=(',',':')))
 records.append({'id':id,'vertices':n,'quadratureFaces':len(quadrature['faces'])//3 if quadrature else None,'sourceClosed':bool(mesh.is_watertight),'sourceOriented':bool(mesh.is_winding_consistent),'sourceVolume':float(mesh.volume),'maxSmoothedOffsetMm':float(np.linalg.norm(D,axis=1).max()/B['definition']['source']['worldUnitsPerMm'])});print(json.dumps(records[-1]),flush=True)
 (out/'report.json').write_text(json.dumps({'method':__doc__,'meshes':records},indent=2)+'\n')
