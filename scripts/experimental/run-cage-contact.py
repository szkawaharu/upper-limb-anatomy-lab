"""Offline sampled actual-bone contact, on a completed material equilibrium.
Sequential tangent planes and penalty continuation approximate unilateral hard
contact. The final test uses the actual closed bones, not the tangent planes.
No source pin is moved. Baseline atlas overlap is allowed point by point.
"""
import sys,json,time,hashlib,shutil
from pathlib import Path
import numpy as np,trimesh
from scipy.spatial.transform import Rotation
from volumetric_reference import Solid
from cage_solver import minimize
from contact_potential import PlaneContact,AugmentedSolid
from bone_sets import model_id,bone_set,DEFINITIONS

file=Path(sys.argv[1]);data=np.load(file);previous=json.loads(file.with_suffix('.json').read_text());p=file.parent;cage=np.load(p/'cage.npz');s=np.load(p/'conforming.npz');assert previous['progress']==1 and previous['allStagesConverged'];assert hashlib.sha256(file.read_bytes()).hexdigest()==previous['snapshotSha256']
identity=model_id(file,previous);coverage=bone_set(identity)
out=p/'contact';out.mkdir(exist_ok=True);output=out/file.stem;clock=time.time();scale=float(cage['worldUnitsPerMm']);center=cage['centerWorld'];rest=cage['nodes'];faces=s['faces'];source=s['points'];ids=s['embeddingVertices'];bary=s['embeddingBarycentric'];fixed=np.flatnonzero(s['classes']>=0);solid=Solid(rest,cage['tets'],fixed,element_weights=cage['weights']);x=data['nodes'].copy();initial_fixed=x[fixed].copy();target=previous['target'];RA=Rotation.from_quat(target['a']['q']).as_matrix();ta=np.array(target['a']['p'])-RA@np.array(target['a']['rest']);surface=lambda nodes:np.einsum('ni,nij->nj',bary,nodes[ids])
code_names=['run-cage-contact.py','contact_potential.py','bone_sets.py','contact-bone-sets.json','cage_solver.py','projected_solid_metric.py','constitutive_metric.py','volumetric_reference.py'];code_hashes={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()for name in code_names};snapshot=out/'code-snapshot';snapshot.mkdir(exist_ok=True)
for name in code_names:shutil.copy2(Path(__file__).with_name(name),snapshot/name)
calfile=Path('docs/validation/v076/source-bones/browser-numeric/rig-calibration.json');cal=json.loads(calfile.read_text());assert scale==cal['scale'];P=np.array([[-1,0,0],[0,0,1],[0,1,0]]);offset=np.array(cal['landmarks']['SC'])/scale-P@np.array([-14.91024,-129.66467,1341.37453]);bones=[]
for entry in coverage['bones']:
 key,name=entry['key'],entry['file']
 path=Path('/tmp/upper-limb-atlas')/name;mesh=trimesh.load_mesh(path,process=True);mesh.vertices=np.asarray(mesh.vertices)@P.T+offset-center/scale;assert mesh.is_watertight
 RB=Rotation.from_quat(target[key]['q']).as_matrix();tb=np.array(target[key]['p'])-RB@np.array(target[key]['rest']);R=RA.T@RB;t=((center@RA.T+ta-tb)@RB-center)/scale
 bones.append({'key':key,'name':name,'mesh':mesh,'R':R,'t':t,'path':path})

def distances(mesh,points):
 inside=mesh.contains(points);use=np.flatnonzero(inside);depth=np.zeros(len(points));normal=np.zeros_like(points)
 for start in range(0,len(use),256):
  idx=use[start:start+256];near,d,face=trimesh.proximity.closest_point(mesh,points[idx]);depth[idx]=d;normal[idx]=(near-points[idx])/np.maximum(d[:,None],1e-30);zero=d<1e-12;normal[idx[zero]]=mesh.face_normals[face[zero]]
 return depth,normal,inside

for b in bones:b['baseline'],_,_=distances(b['mesh'],source)
tri=source[faces];area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2;area_weights=np.zeros(len(source))
for k in range(3):np.add.at(area_weights,faces[:,k],area/3)
area_weights/=area_weights.sum();edges=np.unique(np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0);length=np.linalg.norm(source[edges[:,1]]-source[edges[:,0]],axis=1);penalty=1/float(np.median(length))**2;contact_tolerance_mm=.1

def scan(nodes,make_planes=False):
 points=surface(nodes);rows=[];planes=[]
 for b in bones:
  depth,normal,inside=distances(b['mesh'],points@b['R']+b['t']);increase=depth-b['baseline'];worst=int(increase.argmax());rows.append({'bone':b['name'],'insidePoints':int(inside.sum()),'newInsidePoints':int((inside&(b['baseline']==0)).sum()),'maxBaselineDepthMm':float(b['baseline'].max()),'maxDepthMm':float(depth.max()),'maxDepthIncreaseMm':float(increase.max()),'worstSurfaceVertex':worst,'increaseOverOneMm':int((increase>1).sum()),'increaseOverPointOneMm':int((increase>.1).sum())})
  if make_planes:
   use=np.flatnonzero(inside);n=normal[use]@b['R'].T;boundary=np.einsum('ni,ni->n',n,points[use])+increase[use];no_free=np.all((s['classes'][ids[use]]>=0)|(bary[use]<=1e-9),axis=1)
   conflict=no_free&(increase[use]>1e-5)
   if conflict.any():
    bad=use[conflict];detail={'experimental':True,'accepted':False,'model':identity,'bone':b['name'],'input':str(file),'inputSha256':hashlib.sha256(file.read_bytes()).hexdigest(),'reason':'Fixed boundary contradicts actual bone contact; no pin is moved or removed','vertices':bad.tolist(),'maxDepthIncreaseMm':float(increase[bad].max()),'boneCoverage':coverage};(out/'fixed-boundary-conflict.json').write_text(json.dumps(detail,indent=2)+'\n');raise ValueError(detail['reason'])
   planes.append((use,n,boundary))
 return rows,planes

before,_=scan(x);history=[];reached=False;frame_error=None;independent=file.with_name(file.stem+'-bones.json')
if independent.exists():
 expected=json.loads(independent.read_text());frame_error=max(abs(actual[key]-next(b for b in expected['bones']if b['bone']==actual['bone'])[other])for actual in before for key,other in [('maxBaselineDepthMm','maxRestDepthMm'),('maxDepthMm','maxPosedDepthMm'),('maxDepthIncreaseMm','maxDepthIncreaseMm')]);assert frame_error<1e-7,'Bone/rest/material frame mismatch with independent world-space audit'
for outer in range(8):
 rows,planes=scan(x,True);use=np.concatenate([a[0]for a in planes]);n=np.concatenate([a[1]for a in planes]);boundary=np.concatenate([a[2]for a in planes]);potential=PlaneContact(len(x),ids[use],bary[use],n,boundary,area_weights[use],penalty);augmented=AugmentedSolid(solid,[potential]);trial,report=minimize(augmented,x,bulk=previous['bulk'],iterations=150)
 assert np.array_equal(trial[fixed],initial_fixed);x=trial;rows,_=scan(x);maximum=max(r['maxDepthIncreaseMm']for r in rows);history.append({'outer':outer+1,'penaltyPerMm2':penalty,'planeSamples':len(use),'seconds':time.time()-clock,'actualBones':rows,'minimization':report});print(json.dumps({k:v for k,v in history[-1].items()if k!='minimization'}|{'converged':report['converged'],'residualMaxMm':report['residualMaxMm']}),flush=True)
 if report['converged']and maximum<=contact_tolerance_mm:reached=True;break
 if not report['converged']:break
 penalty*=10

now=surface(x);E=lambda points:np.stack([points[faces[:,k]]-points[faces[:,0]]for k in [1,2]],axis=2);R=E(source);N=E(now);normal0=np.cross(R[:,:,0],R[:,:,1]);normal1=np.cross(N[:,:,0],N[:,:,1]);F=(solid.edges(x)@solid.inverse)[s['elements']];expected=np.linalg.det(F)[:,None]*np.einsum('tij,tj->ti',np.linalg.inv(F).transpose(0,2,1),normal0);vol=lambda y:float(np.einsum('ti,ti->t',y[faces[:,0]],np.cross(y[faces[:,1]],y[faces[:,2]])).sum()/6);er=np.linalg.norm(now[edges[:,1]]-now[edges[:,0]],axis=1)/length;ar=np.linalg.norm(normal1,axis=1)/np.linalg.norm(normal0,axis=1)
result={**previous,'experimental':True,'accepted':False,'contactReached':reached,'contactToleranceMm':contact_tolerance_mm,'contactMethod':'sequential actual closed-bone tangent planes; surface dual-area penalty; final per-material-vertex exact bone depth increase audit','contactBefore':before,'contactHistory':history,'frameConsistencyMaxMm':frame_error,'allStagesConverged':previous['allStagesConverged']and all(h['minimization']['converged']for h in history),'boneContactEvaluated':True,'selfCollisionEvaluated':False,'sourceOpposedFaces':int((np.einsum('ti,ti->t',normal1,expected)<0).sum()),'surfaceVolumeRatio':vol(now)/vol(source),'surfaceEdgeRatio':[float(er.min()),float(er.max())],'surfaceAreaRatio':[float(ar.min()),float(ar.max())],'attachmentMaxMm':float(np.linalg.norm(x[fixed]-initial_fixed,axis=1).max(initial=0)),'seconds':time.time()-clock,'activation':None,'force':None,'contactInputs':{str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in [file,file.with_suffix('.json'),calfile,*[b['path']for b in bones]]},'contactCodeSha256':code_hashes}
result['contactModelId']=identity;result['contactBoneCoverage']=coverage;result['boneDefinitionsSha256']=hashlib.sha256(DEFINITIONS.read_bytes()).hexdigest()
temporary=output.with_suffix('.tmp.npz');np.savez_compressed(temporary,**{**{k:data[k]for k in data.files},'nodes':x,'sourceWorld':(now*scale+center)@RA.T+ta});result['snapshotSha256']=hashlib.sha256(temporary.read_bytes()).hexdigest();temporary.replace(output.with_suffix('.npz'));temporary=output.with_suffix('.tmp.json');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(output.with_suffix('.json'))
print('FINISHED',json.dumps({'contactReached':reached,'maxDepthIncreaseMm':maximum,'surfaceVolumeRatio':result['surfaceVolumeRatio'],'sourceOpposedFaces':result['sourceOpposedFaces'],'seconds':time.time()-clock}),flush=True)
if not reached:raise RuntimeError('Sampled bone contact target not reached; unaccepted result saved')
