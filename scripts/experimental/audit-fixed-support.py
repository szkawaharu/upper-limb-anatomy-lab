"""Detect tissue fixed by interpolation supports beyond its declared mask.
No pose data or deformation score is used to reclassify or remove a footprint.
"""
import sys,json,hashlib
from pathlib import Path
import numpy as np
p=Path(sys.argv[1]);c=np.load(p/'cage.npz');s=np.load(p/'conforming.npz');ids=s['embeddingVertices'];bary=s['embeddingBarycentric'];classes=s['classes'][ids];active=bary>1e-9;unlabelled=~np.isin(s['weights'],[0,1]);all_fixed=np.all((classes>=0)|~active,axis=1);rows=[]
for end in [0,1]:
 mask=unlabelled&np.all((classes==end)|~active,axis=1);which=np.flatnonzero(mask);rows.append({'end':end,'undeclaredRigidVertices':len(which),'firstVertices':which[:20].tolist()})
examples=[]
for i in [12473,28359]:
 if i>=len(ids):continue
 examples.append({'vertex':i,'declaredWeight':float(s['weights'][i]),'controlVertices':ids[i].tolist(),'barycentric':bary[i].tolist(),'supportClasses':classes[i].tolist(),'fixedInfluence':float(bary[i][classes[i]>=0].sum())})
report={'experimental':True,'accepted':False,'sourceDirectory':str(p),'inputSha256':{name:hashlib.sha256((p/name).read_bytes()).hexdigest()for name in ['cage.npz','conforming.npz']},'surfaceVertices':len(ids),'fixedControls':int((s['classes']>=0).sum()),'declaredFixedSurfaceVertices':int((~unlabelled).sum()),'unlabelledWithNoFreeSupport':int((all_fixed&unlabelled).sum()),'byBone':rows,'examples':examples,'classificationTolerance':'1e-9 barycentric roundoff screening; not an anatomical or distance acceptance threshold'}
(p/'fixed-support-audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
