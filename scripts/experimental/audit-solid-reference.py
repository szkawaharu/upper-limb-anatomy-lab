"""Independent geometry audit of a saved offline source-solid pose.
Edge/face/volume measurements do not certify collision freedom or anatomy.
"""
import sys,json,hashlib
from pathlib import Path
import numpy as np
f=Path(sys.argv[1]);out=Path(sys.argv[2]);a=np.load(f);r=a['referenceWorld'];p=a['world'];t=a['tets'];faces=a['sourceFaces'].reshape(-1,3);scale=.00555

def signed_volumes(points):
 q=points[t];return np.einsum('ti,ti->t',q[:,1]-q[:,0],np.cross(q[:,2]-q[:,0],q[:,3]-q[:,0]))/6
restv=signed_volumes(r);nowv=signed_volumes(p);J=nowv/restv
assert np.all(restv>0)
ids=np.unique(faces);assert np.all(np.isfinite(p))
edge=np.unique(np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0)
le=lambda points:np.linalg.norm(points[edge[:,1]]-points[edge[:,0]],axis=1)
area=lambda points:np.linalg.norm(np.cross(points[faces[:,1]]-points[faces[:,0]],points[faces[:,2]]-points[faces[:,0]]),axis=1)/2
edge_ratio=le(p)/le(r);area_ratio=area(p)/area(r)
# Float32 reference AND deformed points reflect actual WebGL input precision.
rv32=signed_volumes(r.astype(np.float32).astype(float));pv32=signed_volumes(p.astype(np.float32).astype(float));J32=pv32/rv32
boundary_edges=np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1);_,multiplicity=np.unique(boundary_edges,axis=0,return_counts=True)
surface_volume=lambda points:float(np.einsum('ti,ti->t',points[faces[:,0]]-points[ids].mean(0),np.cross(points[faces[:,1]]-points[ids].mean(0),points[faces[:,2]]-points[ids].mean(0))).sum()/6)
report={'experimental':True,'accepted':False,'input':str(f),'inputSha256':hashlib.sha256(f.read_bytes()).hexdigest(),'worldUnitsPerMm':scale,'sourceSurfaceVertices':len(ids),'sourceSurfaceTriangles':len(faces),'sourceBoundaryEdgeMultiplicities':{str(k):int((multiplicity==k).sum())for k in np.unique(multiplicity)},'tetrahedra':len(t),'nonfiniteCoordinates':int((~np.isfinite(p)).sum()),'nonPositiveElements':int((J<=0).sum()),'JQuantiles':np.quantile(J,[0,.01,.5,.99,1]).tolist(),'tetVolumeRatio':float(nowv.sum()/restv.sum()),'surfaceVolumeRatio':surface_volume(p)/surface_volume(r),'tetVsSurfaceDeformedVolumeRelativeError':float(abs(nowv.sum()-surface_volume(p))/abs(nowv.sum())),'edgeRatioMinMax':[float(edge_ratio.min()),float(edge_ratio.max())],'areaRatioMinMax':[float(area_ratio.min()),float(area_ratio.max())],'float32':{'referenceNonPositiveElements':int((rv32<=0).sum()),'deformedNonPositiveElements':int((pv32<=0).sum()),'JQuantiles':np.quantile(J32,[0,.01,.5,.99,1]).tolist()},'worstElements':[],'boneCollisionEvaluated':False,'selfCollisionEvaluated':False,'medicalValidation':False}
for k in [J.argmin(),J.argmax()]:
 q=r[t[k]];E=np.column_stack([q[i]-q[0]for i in [1,2,3]])
 report['worstElements'].append({'tetrahedron':int(k),'vertices':t[k].tolist(),'J':float(J[k]),'referenceVolumeMm3':float(restv[k]/scale**3),'referenceCondition':float(np.linalg.cond(E)),'fixedVertices':np.isin(t[k],a['fixed']).tolist()})
out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
