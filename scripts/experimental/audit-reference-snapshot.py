"""Recompute snapshot metrics; never inherit them from its initial report.
This audits saved geometry without rerunning or claiming a new equilibrium.
"""
import sys,json,hashlib
from pathlib import Path
import numpy as np
from material_pose import boundary_targets,material_to_world
root,file,out=map(Path,sys.argv[1:4]);r=json.loads(file.with_suffix('.json').read_text());assert r['progress']==1 and r['allStagesConverged']and r['contactReached'];assert hashlib.sha256(file.read_bytes()).hexdigest()==r['snapshotSha256']
a=np.load(root/'cage.npz');s=np.load(root/'conforming.npz');data=np.load(file);rest=a['nodes'];x=data['nodes'];T=a['tets'];FACES=s['faces'];assert np.isfinite(x).all();scale=float(a['worldUnitsPerMm']);center=a['centerWorld']
edge=lambda v:np.stack([v[T[:,k]]-v[T[:,0]]for k in [1,2,3]],axis=2)
F=edge(x)@np.linalg.inv(edge(rest));J=np.linalg.det(F);source=s['points'];points=np.einsum('ni,nij->nj',s['embeddingBarycentric'],x[s['embeddingVertices']]);expected_world=material_to_world(points,r['target'],center,scale);world_error=float(np.linalg.norm(expected_world-data['sourceWorld'],axis=1).max()/scale);assert world_error<1e-7
volume=lambda p:float(np.einsum('ti,ti->t',p[FACES[:,0]],np.cross(p[FACES[:,1]],p[FACES[:,2]])).sum()/6)
ratio=volume(points)/volume(source);assert abs(ratio-r['surfaceVolumeRatio'])<1e-8
fixed=np.flatnonzero(s['classes']>=0);prescribed=boundary_targets(rest,fixed,s['classes'],r['target'],center,scale);pin_error=float(np.linalg.norm(x[fixed]-prescribed,axis=1).max(initial=0))
edges=np.unique(np.sort(np.concatenate([FACES[:,[0,1]],FACES[:,[1,2]],FACES[:,[2,0]]]),axis=1),axis=0);er=np.linalg.norm(points[edges[:,1]]-points[edges[:,0]],axis=1)/np.linalg.norm(source[edges[:,1]]-source[edges[:,0]],axis=1)
n0=np.cross(source[FACES[:,1]]-source[FACES[:,0]],source[FACES[:,2]]-source[FACES[:,0]]);n=np.cross(points[FACES[:,1]]-points[FACES[:,0]],points[FACES[:,2]]-points[FACES[:,0]]);area=np.linalg.norm(n,axis=1)/np.linalg.norm(n0,axis=1);selected=F[s['elements']];expected=np.linalg.det(selected)[:,None]*np.einsum('tij,tj->ti',np.linalg.inv(selected).transpose(0,2,1),n0);cells=np.unique(s['elements'])
physical=a['sourceMaterialVolumes']if 'sourceMaterialVolumes'in a else None
physical_metrics={'positiveSourceIntersectionJRange':None,'mappedPhysicalVolumeRatio':None,'surfaceIntegralDifference':None}
if physical is not None:
 mask=physical>0;physical_ratio=float(np.dot(physical,J)/physical.sum());physical_metrics={'positiveSourceIntersectionJRange':[float(J[mask].min()),float(J[mask].max())],'positiveSourceIntersectionCells':int(mask.sum()),'mappedPhysicalVolumeRatio':physical_ratio,'surfaceIntegralDifference':abs(physical_ratio-ratio),'materialCellMask':'strictly positive source/cage intersection volume; no size threshold; floating geometric classification'}
assert np.isfinite(J).all()and np.isfinite(er).all()and np.isfinite(area).all()
result={**physical_metrics,'experimental':True,'accepted':False,'postAuditWithoutResolving':True,'modelId':r['contactModelId'],'target':r['target'],'bulk':r['bulk'],'surfaceSignedVolumeRatio':ratio,'attachmentMaxMm':pin_error,'worldFrameMaxErrorMm':world_error,'cageJRange':[float(J.min()),float(J.max())],'surfaceAdjacentJRange':[float(J[cells].min()),float(J[cells].max())],'cageNonPositive':int((J<=0).sum()),'surfaceEdgeRatio':[float(er.min()),float(er.max())],'surfaceAreaRatio':[float(area.min()),float(area.max())],'sourceOpposedFaces':int((np.einsum('ti,ti->t',n,expected)<0).sum()),'actualDistanceResidualMm':r['actualDistanceResidualMm'],'lastSolvedActualBoneDepth':r['history'][-1]['bones'],'rawSnapshot':str(file),'rawSnapshotSha256':r['snapshotSha256'],'rawReportSha256':hashlib.sha256(file.with_suffix('.json').read_bytes()).hexdigest(),'inputsSha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in [root/'cage.npz',root/'conforming.npz']},'codeSha256':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()for name in ['audit-reference-snapshot.py','material_pose.py']},'limitations':['Metrics are recomputed from saved vertices; inherited initial-report geometry metrics are not used','Signed surface volume is not anatomical/global self-collision certification','Bone depth is from the nonlinear solve; independent world-space bone audit is separate','Not a production GPU or mobile test']}
out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in ['target','inputsSha256','limitations','codeSha256','lastSolvedActualBoneDepth']}))
