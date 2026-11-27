"""Split source triangles along embedding tetrahedra without moving the source
surface. Each output triangle lies in one affine cell. Original vertices are
retained; source-face and barycentric provenance is recorded for new vertices.
"""
import itertools
import numpy as np


def conform(cage,points,faces,weights=None,callback=None):
 points=np.asarray(points);faces=np.asarray(faces).reshape(-1,3);original=len(points);vertices=points.tolist();origins=[(int(i),-1,[1.,0.,0.])for i in range(original)];new_weights=list(weights)if weights is not None else None;lookup={};triangles=[];source_faces=[];cells=[];max_merge=0.
 _,_,source_cells=cage.embed(points)
 inv=np.linalg.inv(np.stack([cage.nodes[cage.tets[:,k]]-cage.nodes[cage.tets[:,0]]for k in [1,2,3]],axis=2));anchors=cage.nodes[cage.tets[:,0]]
 def vertex(face_id,face,u):
  nonlocal max_merge
  corner=int(np.argmax(u))
  if abs(u[corner]-1)<1e-10:return int(face[corner])
  p=u@points[face];zero=np.flatnonzero(np.abs(u)<1e-10)
  if len(zero):
   edge=[k for k in range(3)if k!=zero[0]];a,b=sorted(int(face[k])for k in edge);ratio=float(u[np.flatnonzero(face==b)[0]]);key=('edge',a,b,round(ratio,10))
  else:key=('face',face_id,*np.round(u,10))
  if key in lookup:
   i=lookup[key];error=float(np.linalg.norm(np.asarray(vertices[i])-p));max_merge=max(max_merge,error)
   if error>1e-7:raise ValueError('Inconsistent conforming shared vertex')
   return i
  i=len(vertices);lookup[key]=i;vertices.append(p.tolist());origins.append((-1,face_id,u.tolist()))
  if weights is not None:
   w=np.asarray(weights)[face];new_weights.append(float(w[0])if np.all(w==w[0])else float(u@w))
  return i
 for face_id,face in enumerate(faces):
  src=points[face];same=source_cells[face]
  if np.all(same==same[0]):triangles.append(face.tolist());source_faces.append(face_id);cells.append(int(same[0]));continue
  low=np.floor((src.min(0)-cage.origin)/cage.h).astype(int);high=np.floor((src.max(0)-cage.origin)/cage.h).astype(int)
  for cell in itertools.product(*[range(low[k],high[k]+1)for k in range(3)]):
   first=int(np.ravel_multi_index(cell,cage.cells)*6)
   for tet in range(first,first+6):
    tail=(src-anchors[tet])@inv[tet].T;coordinates=np.column_stack([1-tail.sum(1),tail]);poly=[np.eye(3)[k]for k in range(3)]
    for plane in range(4):
     if not poly:break
     clipped=[]
     for j,b in enumerate(poly):
      a=poly[j-1];sa=float(a@coordinates[:,plane]);sb=float(b@coordinates[:,plane]);ina=sa>=-1e-12;inb=sb>=-1e-12
      if ina!=inb:clipped.append(a+(b-a)*(sa/(sa-sb)))
      if inb:clipped.append(b)
     poly=clipped
    if len(poly)<3:continue
    for j in range(1,len(poly)-1):
     u=np.array([poly[0],poly[j],poly[j+1]]);area=np.linalg.det(np.column_stack([u[1,1:]-u[0,1:],u[2,1:]-u[0,1:]]))
     if abs(area)<1e-12:continue
     tri=[vertex(face_id,face,x)for x in u]
     if len(set(tri))<3:continue
     triangles.append(tri);source_faces.append(face_id);cells.append(tet)
  if callback and face_id%2000==0:callback(face_id)
 v=np.array(vertices);f=np.array(triangles,np.int32);p=v[f];area=np.linalg.norm(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]),axis=1);orig_p=points[faces];original_area=np.linalg.norm(np.cross(orig_p[:,1]-orig_p[:,0],orig_p[:,2]-orig_p[:,0]),axis=1);area_sum=np.bincount(source_faces,weights=area,minlength=len(faces));relative=np.abs(area_sum-original_area)/np.maximum(original_area,1e-30)
 if relative.max()>1e-6:raise ValueError('Source triangle area coverage changed: '+str(relative.max()))
 assert np.array_equal(v[:original],points),'Original vertices changed'
 return {'points':v,'faces':f,'weights':np.array(new_weights)if weights is not None else None,'sourceTriangle':np.array(source_faces),'element':np.array(cells),'newVertexProvenance':origins[original:],'sourceAreaMaxRelativeError':float(relative.max()),'maxMergeDistanceMm':max_merge}
