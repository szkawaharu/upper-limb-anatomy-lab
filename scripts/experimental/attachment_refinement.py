"""Unlock artificial fixed edges/faces, preserving declared attachment facets.
P1 interpolation otherwise makes a face rigid whenever its three controls are
fixed, even if that triangle was never a declared tissue attachment. Neutral
centroid refinement adds free controls on those artificial simplices. It is a
discretization change, independent of pose, contact scores or muscle identity.
"""
import itertools
import numpy as np

def refine(nodes,tets,classes,fixed_facets,fixed_edges,needed=None):
 nodes=np.asarray(nodes);tets=np.asarray(tets);classes=np.asarray(classes);points=nodes.tolist();labels=classes.tolist();declared_faces={tuple(sorted(f))for f in fixed_facets};declared_edges={tuple(sorted(e))for e in fixed_edges};edges={};faces={};tet_faces=[]
 for t in tets:
  for e in itertools.combinations(t,2):
   key=tuple(sorted(e))
   if key not in declared_edges and classes[e[0]]>=0 and classes[e[0]]==classes[e[1]]and(needed is None or key in needed['edges']):edges.setdefault(key,None)
  tet_faces.append([tuple(sorted(f))for f in itertools.combinations(t,3)])
 for key in edges:edges[key]=len(points);points.append(nodes[list(key)].mean(0).tolist());labels.append(-1)
 for fs in tet_faces:
  for key in fs:
   if key in faces:continue
   bad_edge=any(tuple(sorted(e))in edges for e in itertools.combinations(key,2));same=classes[key[0]]>=0 and np.all(classes[list(key)]==classes[key[0]])
   if key not in declared_faces and (bad_edge or (same and(needed is None or key in needed['faces']))):faces[key]=None
 for key in faces:faces[key]=len(points);points.append(nodes[list(key)].mean(0).tolist());labels.append(-1)
 cells=[];parents=[]
 for parent,(t,fs)in enumerate(zip(tets,tet_faces)):
  if not any(f in faces for f in fs)and not(needed is not None and parent in needed['cells']):cells.append(t.tolist());parents.append(parent);continue
  center=len(points);points.append(nodes[t].mean(0).tolist());labels.append(-1)
  for f in fs:
   if f not in faces:cells.append([center,*f]);parents.append(parent);continue
   polygon=[]
   for k in range(3):
    a,b=f[k],f[(k+1)%3];polygon.append(a);key=tuple(sorted([a,b]))
    if key in edges:polygon.append(edges[key])
   for k in range(len(polygon)):cells.append([center,faces[f],polygon[k],polygon[(k+1)%len(polygon)]]);parents.append(parent)
 p=np.array(points);c=np.array(cells,np.int32);E=lambda p,t:np.stack([p[t[:,k]]-p[t[:,0]]for k in [1,2,3]],axis=2);det=np.linalg.det(E(p,c));assert np.all(det!=0)
 reverse=det<0;c[reverse,1],c[reverse,2]=c[reverse,2].copy(),c[reverse,1].copy();old_volume=float(np.linalg.det(E(nodes,tets)).sum()/6);new_volume=float(np.abs(det).sum()/6);assert abs(new_volume-old_volume)<1e-10*max(1.,old_volume)
 return p,c,np.array(labels,np.int8),{'addedFreeEdgeControls':len(edges),'addedFreeFaceControls':len(faces),'addedFreeCellControls':len(points)-len(nodes)-len(edges)-len(faces),'referenceVolumeDifference':new_volume-old_volume,'poseIndependent':True,'originalControlPointsChanged':False},np.array(parents)
