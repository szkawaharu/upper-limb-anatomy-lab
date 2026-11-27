"""Independently intersect a cage cell with the original fine material mesh.
No integration weights are changed by this diagnostic. An unresolved sliver
is reported, not silently counted as a certified zero-volume intersection.
"""
import sys,json,hashlib,time
from pathlib import Path
import numpy as np
from rtree import index
from tet_intersection import intersection

root,material_path,out=map(Path,sys.argv[1:4])
cell_ids=json.loads(sys.argv[4])
a=np.load(root/'cage.npz');m=np.load(material_path)
scale=float(a['worldUnitsPerMm']);center=a['centerWorld']
fine=(m['nodes']-center)/scale;tets=m['tets'];tri=fine[tets]
low,high=tri.min(1),tri.max(1)
prop=index.Property();prop.dimension=3
tree=index.Index(((i,tuple(np.r_[low[i],high[i]]),None)for i in range(len(tets))),properties=prop)
source_volume=float(np.linalg.det(np.stack([tri[:,k]-tri[:,0]for k in [1,2,3]],axis=2)).sum()/6)
rows=[];clock=time.time()
for cell in cell_ids:
 tet=a['nodes'][a['tets'][cell]];bounds=np.r_[tet.min(0),tet.max(0)];candidates=list(tree.intersection(tuple(bounds)))
 measured=0.;statuses={};unresolved=[]
 for i in candidates:
  value,report=intersection(tet,tri[i]);measured+=value;status=report['status'];statuses[status]=statuses.get(status,0)+1
  if status=='rank-unresolved-or-contact':unresolved.append(i)
 cage_volume=float(abs(np.linalg.det(np.stack([tet[k]-tet[0]for k in [1,2,3]],axis=1)))/6)
 # Existing positive regularization and normalization remain included; this
 # is an equivalent material-volume contribution, not a measured density.
 equivalent=float(a['weights'][cell]*source_volume)
 row={'cell':cell,'cageVolumeMm3':cage_volume,'fineTetrahedraCandidates':len(candidates),'intersectionVolumeMm3':measured,'existingWeightEquivalentVolumeMm3':equivalent,'existingOverIntersectionRatio':equivalent/measured if measured>0 else None,'statuses':statuses,'rankUnresolvedOrContactFineTets':unresolved,'seconds':time.time()-clock}
 rows.append(row);print(json.dumps(row),flush=True)
result={'experimental':True,'accepted':False,'weightsChanged':False,'codeSha256':{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()for name in ['audit-cell-integration.py','tet_intersection.py']},'method':'convex fine-tet/cage-tet intersection; not exact rational arithmetic','sourceVolumeMm3':source_volume,'rows':rows,'inputSha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in [root/'cage.npz',material_path]}}
out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n')
