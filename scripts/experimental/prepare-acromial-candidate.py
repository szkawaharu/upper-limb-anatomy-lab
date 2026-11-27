"""Reclassify the entire legacy origin from neutral anatomical regions only.
This is an explicit unadopted boundary-data experiment, not per-pose pin removal.
The deltoid source mesh and every vertex survive unchanged. A broad near-bone
region is screened from the lateral acromial wall; fascia is not relabelled as
measured bone enthesis. Stop-side estimated tissue mask remains unverified.
"""
import sys,json,hashlib
from pathlib import Path
import numpy as np,trimesh
from scipy.spatial import cKDTree

source=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
definition_file=Path(__file__).with_name('acromial-region-candidate.json');D=json.loads(definition_file.read_text());calfile=Path('docs/validation/v076/source-bones/browser-numeric/rig-calibration.json');cal=json.loads(calfile.read_text());scale=cal['scale'];L={k:np.array(v)/scale for k,v in cal['landmarks'].items()}
B=next(b for b in json.loads(source.read_text())if b['spec']['id']==D['id']);points=np.array(B['rest']).reshape(-1,3)/scale;faces=np.array(B['faces']).reshape(-1,3);old=np.array(B['weights']);edges=np.unique(np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0);spacing=float(np.median(np.linalg.norm(points[edges[:,0]]-points[edges[:,1]],axis=1)))
bonefile=Path('/tmp/upper-limb-atlas/FMA13395.stl');bone=trimesh.load_mesh(bonefile,process=True);P=np.array([[-1,0,0],[0,0,1],[0,1,0]]);offset=L['SC']-P@np.array([-14.91024,-129.66467,1341.37453]);bone.vertices=np.array(bone.vertices)@P.T+offset
dist=lambda names:np.linalg.norm(np.asarray(bone.vertices)[:,None,:]-np.array([L[k]for k in names])[None,:,:],axis=2).min(1)
region=dist(D['includedLandmarks'])<dist(D['excludedLandmarks']);outward=L[D['outwardAxis'][1]]-L[D['outwardAxis'][0]];outward/=np.linalg.norm(outward);region&=np.asarray(bone.vertex_normals)@outward>0
roi=trimesh.Trimesh(vertices=bone.vertices.copy(),faces=bone.faces[region[bone.faces].all(1)],process=False);roi.remove_unreferenced_vertices();assert len(roi.faces)>0
axis=L[D['proximalAxis'][1]]-L[D['proximalAxis'][0]];axis/=np.linalg.norm(axis);allowed=np.flatnonzero(((points-L[D['proximalHalfSpaceOrigin']])@axis>=0)&(old!=1));distance=np.full(len(points),np.inf);nearest=np.zeros_like(points)
for i in range(0,len(allowed),256):
 ids=allowed[i:i+256];q,d,_=trimesh.proximity.closest_point(roi,points[ids]);distance[ids]=d;nearest[ids]=q
minimum=float(distance.min());resolution=spacing+float(np.median(bone.edges_unique_length));origin=np.flatnonzero(distance<=minimum+resolution);assert len(origin)>=3
weights=np.clip(old,2**-24,1-2**-24);weights[old==1]=1;weights[origin]=0
assert not np.any((old==1)&(weights==0))
report={'experimental':True,'adopted':False,'definition':D,'inputs':{str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in [source,definition_file,calfile,bonefile]},'restOnly':True,'legacyOriginVertices':int((old==0).sum()),'candidateOriginVertices':len(origin),'retainedLegacyOriginVertices':int(((weights==0)&(old==0)).sum()),'legacyOriginReclassified':int(((weights!=0)&(old==0)).sum()),'candidateIndices':origin.tolist(),'originGapMm':[float(distance[origin].min()),float(distance[origin].max())],'resolutionBandMm':resolution,'muscleMedianEdgeMm':spacing,'originCentroidAxisMm':points[origin].mean(0).tolist(),'boneROIVertices':len(roi.vertices),'boneROITriangles':len(roi.faces),'oppositeFootprintUnchanged':True,'verticesAndFacesChanged':False,'notAClaimOfClinicalAttachmentAccuracy':True}
B['weights']=weights.tolist();B['material']=weights.tolist();B['boundaryCandidate']=report;(out/'bindings.json').write_text(json.dumps([B],separators=(',',':')))
np.savez_compressed(out/'boundary.npz',pointsWorld=points*scale,faces=faces,weights=weights,legacyWeights=old,origin=origin,boneROIPointsWorld=roi.vertices*scale,boneROIFaces=roi.faces)
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in ['inputs','candidateIndices','definition']}))
