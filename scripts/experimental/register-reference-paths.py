"""Cross-atlas path diagnostic, not a new muscle footprint definition.
Scapula: TS/AA/AI proper similarity. Humerus: head/elbow length plus transverse
direction. The fitted condylar centres are not treated as epicondylar width.
"""
import sys,json,hashlib
from pathlib import Path
import numpy as np,trimesh
from landmark_registration import similarity,bone_similarity,apply

src=Path('docs/validation/v076/reference-model/subject-paths.json');calfile=Path('docs/validation/v076/source-bones/browser-numeric/rig-calibration.json')
data=json.loads(src.read_text());cal=json.loads(calfile.read_text());L={k:np.array(v)/cal['scale']for k,v in cal['landmarks'].items()};markers={x['name']:np.array(x['locationMetres'])*1000 for x in data['markers']}
P=np.array([[-1,0,0],[0,0,1],[0,1,0]]);offset=L['SC']-P@np.array([-14.91024,-129.66467,1341.37453])
reg={'scapula':similarity(np.array([markers[x]for x in ['ts','aa','ai']]),np.array([L[x]for x in ['TS','AA','AI']])),'humerus':bone_similarity(markers['HumHead_Center'],markers['centelbow'],markers['EpL']-markers['EpM'],L['GH'],L['EL'],L['EL_LAT']-L['EL_MED'])}
meshes={};paths=[]
for name,file in [('scapula','FMA13395.stl'),('humerus','FMA23130.stl')]:
 m=trimesh.load_mesh(Path('/tmp/upper-limb-atlas')/file,process=True);m.vertices=np.array(m.vertices)@P.T+offset;meshes[name]=m
for muscle in data['muscles']:
 row={'name':muscle['name'],'points':[]}
 for p in muscle['points']:
  body=p['bodyFrame'].split('/')[-1]
  if body not in reg:continue
  point=apply(reg[body],np.array(p['locationMetres'])*1000);near,gap,face=trimesh.proximity.closest_point(meshes[body],point[None,:])
  row['points'].append({'name':p['name'],'body':body,'positionAxisMm':point.tolist(),'sourceLocalMetres':p['locationMetres'],'nearestBonePointAxisMm':near[0].tolist(),'nearestBoneFace':int(face[0]),'boneGapMm':float(gap[0])})
 if row['points']:paths.append(row)
serial=lambda r:{k:v.tolist()if isinstance(v,np.ndarray)else v for k,v in r.items()}
report={'experimental':True,'adoptedIntoApplication':False,'definesAttachmentFootprints':False,'units':'render world axes, mm (offset inherited from calibrated atlas)','source':{k:data[k]for k in ['sourceRepository','sourceCommit','sourcePath','modelSha256','license','copyright','attribution']},'inputSha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in [src,calfile]},'registrations':{name:serial(r)for name,r in reg.items()},'limitations':['Different atlas individuals; geometry is not clinically registered','Three estimated scapular landmarks; residuals are not independent holdout validation','Humeral transverse direction from fitted condylar centres; no false epicondylar-width scaling','Muscle line paths cannot independently define broad attachment masks','Paths on unregistered bones have been omitted'],'paths':paths}
out=Path(sys.argv[1]);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'registration':report['registrations'],'deltoid':[p for p in paths if p['name'].startswith('Deltoideus')]}))
