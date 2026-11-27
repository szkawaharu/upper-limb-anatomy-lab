# SPDX-License-Identifier: GPL-3.0-or-later
"""OFFLINE MeshFix candidate with explicit atlas deviation audit.
No component removal/joining is enabled. The immutable original is retained;
this is an unadopted, independently meshed reference, not a pose correction.
The GPL repair library is not linked to the web application.
"""
import sys,json,hashlib,time,importlib.metadata,shutil
from pathlib import Path
import numpy as np,trimesh,pymeshfix
from scipy.spatial import cKDTree

src=Path(sys.argv[1]);id=sys.argv[2];out=Path(sys.argv[3]);out.mkdir(parents=True,exist_ok=True);clock=time.time();B=next(b for b in json.loads(src.read_text())if b['spec']['id']==id);scale=B['definition']['source']['worldUnitsPerMm'];p=np.array(B['rest']).reshape(-1,3)/scale;f=np.array(B['faces']).reshape(-1,3);before=trimesh.Trimesh(p,f,process=False);components_before=len(before.split(only_watertight=False));fix=pymeshfix.MeshFix(p.copy(),f.copy());fix.repair(verbose=False,joincomp=False,remove_smallest_components=False);q=fix.points.copy();g=fix.faces
# Undo only an exact Float32 round trip for retained vertices. Genuine
# changed/new repair points are not snapped, and source tissue is not jiggled.
_,neighbor=cKDTree(p).query(q);round_trip=np.all(q==p[neighbor].astype(np.float32).astype(float),axis=1);q[round_trip]=p[neighbor[round_trip]]
after=trimesh.Trimesh(q,g,process=False);components_after=len(after.split(only_watertight=False));spacing=float(np.median(before.edges_unique_length));gaps=[];nearest_faces=[]
for mesh,points in [(after,p),(before,q)]:
 distance=[];near=[]
 for start in range(0,len(points),256):
  _,d,tri=trimesh.proximity.closest_point(mesh,points[start:start+256]);distance.extend(d);near.extend(tri)
 gaps.append(np.array(distance));nearest_faces.append(np.array(near))
d,nearest=cKDTree(p).query(q);original_exact=d==0;old=np.array(B['weights']);tri=f[nearest_faces[1]];closest,_,_=trimesh.proximity.closest_point(before,q);bary=trimesh.triangles.points_to_barycentric(p[tri],closest);w=np.clip(np.einsum('ni,ni->n',bary,old[tri]),np.nextafter(np.float32(0),np.float32(1)),1-2**-24)
for end in [0,1]:w[np.all(old[tri]==end,axis=1)]=end
w[original_exact]=old[nearest[original_exact]];w=w.astype(np.float32)
report={'experimental':True,'adoptedIntoApplication':False,'id':id,'sourceInput':str(src),'inputSha256':hashlib.sha256(src.read_bytes()).hexdigest(),'repairLibrary':'pymeshfix','version':pymeshfix.__version__,'license':'GPL-3.0-or-later; offline only','joinComponents':False,'removeSmallestComponents':False,'sourceVertices':len(p),'sourceTriangles':len(f),'repairedVertices':len(q),'repairedTriangles':len(g),'componentsBefore':components_before,'componentsAfter':components_after,'closedBefore':before.is_watertight,'closedAfter':after.is_watertight,'orientedAfter':after.is_winding_consistent,'eulerBefore':before.euler_number,'eulerAfter':after.euler_number,'originalExactVerticesRetained':int(original_exact.sum()),'sourceMedianEdgeMm':spacing,'directedVertexSurfaceDistanceQuantilesMm':[np.quantile(a,[0,.5,.99,1]).tolist()for a in gaps],'directedVertexSurfaceDistanceMaxMm':max(float(a.max())for a in gaps),'distanceIsNotContinuousHausdorffCertification':True,'weightTransfer':'exact original vertex where unchanged, otherwise original nearest-triangle barycentric field; original footprint remains an estimate','seconds':time.time()-clock}
report['geometricScreeningPassed']=bool(after.is_watertight and after.is_winding_consistent and components_before==components_after and report['directedVertexSurfaceDistanceMaxMm']<=spacing)
(out/'repair-audit.json').write_text(json.dumps(report,indent=2)+'\n');np.savez_compressed(out/'original.npz',points=p,faces=f,weights=old);np.savez_compressed(out/'candidate.npz',points=q,faces=g,weights=w,sourceNearestVertex=nearest,sourceVertexExact=original_exact,sourceNearestFace=nearest_faces[1],sourceBarycentric=bary)
metadata=importlib.metadata.distribution('pymeshfix')
for path in metadata.files:
 if str(path).upper().endswith('LICENSE'):shutil.copy2(metadata.locate_file(path),out/'repair-library-license.txt');break
if not report['geometricScreeningPassed']:raise RuntimeError('Repair candidate failed atlas deviation/component screening; do not adopt')
definition={'source':{**B['definition']['source'],'repairCandidate':{'method':'MeshFix without component removal or joining','audit':'repair-audit.json','medicalAccuracyVerified':False}}};result={'spec':B['spec'],'rest':(q*scale).ravel().tolist(),'faces':g.ravel().tolist(),'weights':w.tolist(),'sourceVertexIds':np.where(original_exact,np.array(B['sourceVertexIds'])[nearest],-1).tolist(),'definition':definition,'meshPreparation':{**B['meshPreparation'],'closed':after.is_watertight,'oriented':after.is_winding_consistent,'originalGeometryRetained':False,'repairCandidate':report}}
(out/'bindings.json').write_text(json.dumps([result],separators=(',',':')));print(json.dumps(report),flush=True)
