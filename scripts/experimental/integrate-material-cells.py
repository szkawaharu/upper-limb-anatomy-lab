"""Offline source-tet/cage-tet integration, independent of pose or muscle name.
Preserves existing parent-cell background density; replaces four-point source
sampling only. Does not move a source point or change an attachment mask.
"""
import sys,json,time,hashlib,shutil
from pathlib import Path
import numpy as np
from rtree import index
from tet_intersection import prepare,intersections
from constrained_cage import ConstrainedCage

root,material_path,out,parent=map(Path,sys.argv[1:5])
out.mkdir(parents=True,exist_ok=True)
a=np.load(root/'cage.npz');m=np.load(material_path);old=np.load(parent/'cage.npz')
meta=json.loads((root/'preparation.json').read_text());original=json.loads((parent/'preparation.json').read_text())
assert hashlib.sha256((parent/'cage.npz').read_bytes()).hexdigest()==meta['parentInputs']['cage.npz']
scale=float(a['worldUnitsPerMm']);center=a['centerWorld'];fine=(m['nodes']-center)/scale
fine_tets=fine[m['tets']];cage_tets=a['nodes'][a['tets']]
V=lambda t:np.linalg.det(np.stack([t[:,k]-t[:,0]for k in [1,2,3]],axis=2))/6
fine_volume=V(fine_tets);cell_volume=V(cage_tets);assert np.all(fine_volume>0)and np.all(cell_volume>0)
prepared_fine=prepare(fine_tets);prepared_cells=prepare(cage_tets)
code_files=['integrate-material-cells.py','tet_intersection.py'];code_hashes={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()for name in code_files}
(out/'code-snapshot').mkdir(exist_ok=True)
for name in code_files:shutil.copy2(Path(__file__).with_name(name),out/'code-snapshot'/name)
lo,hi=fine_tets.min(1),fine_tets.max(1);prop=index.Property();prop.dimension=3
tree=index.Index(((i,tuple(np.r_[lo[i],hi[i]]),None)for i in range(len(fine_tets))),properties=prop)
physical=np.zeros(len(cage_tets));coverage=np.zeros(len(fine_tets));statuses={};clock=time.time();pairs=0
for i,tet in enumerate(cage_tets):
 candidates=list(tree.intersection(tuple(np.r_[tet.min(0),tet.max(0)])))
 try:volumes,counts=intersections(prepared_cells[i],[prepared_fine[j]for j in candidates])
 except Exception:
  np.savez_compressed(out/'failed-cell.npz',cageCell=i,candidates=candidates,cagePoints=tet,finePoints=fine_tets[candidates]);raise
 physical[i]=volumes.sum();coverage[candidates]+=volumes
 for status,count in counts.items():statuses[status]=statuses.get(status,0)+count
 pairs+=len(candidates)
 if (i+1)%500==0 or i+1==len(cage_tets):
  print(json.dumps({'cells':i+1,'totalCells':len(cage_tets),'pairs':pairs,'seconds':time.time()-clock}),flush=True)
  np.savez_compressed(out/'integration-progress.npz',physicalVolumes=physical,coveredFineVolumes=coverage,cellsCompleted=i+1)
error=coverage-fine_volume;total=float(fine_volume.sum());relative=float(abs(physical.sum()-total)/total);max_error=float(np.max(np.abs(error)));maximum=float(fine_volume.max())
report={'experimental':True,'accepted':False,'method':'convex fine-tet/cage-tet intersections; floating point, not exact rational arithmetic','sourceVolumeMm3':total,'intersectionVolumeMm3':float(physical.sum()),'relativeVolumeCoverageError':relative,'maxFineTetCoverageErrorMm3':max_error,'maxFineTetVolumeMm3':maximum,'sumAbsoluteFineCoverageErrorMm3':float(np.abs(error).sum()),'physicalCells':int((physical>0).sum()),'totalCells':len(cell_volume),'pairs':pairs,'statuses':statuses,'seconds':time.time()-clock,'codeSha256':code_hashes,'sourcePointsChanged':False,'sourceAttachmentMasksChanged':False,'inputsSha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in [root/'cage.npz',root/'conforming.npz',root/'preparation.json',material_path,parent/'cage.npz',parent/'preparation.json']}}
(out/'integration.json').write_text(json.dumps(report,indent=2)+'\n')
# Numerical coverage checks, not a physiological threshold. Keep the raw
# errors and do not emit replacement weights on a coverage failure.
assert relative<1e-8 and max_error<1e-6*maximum, 'Source material coverage mismatch; replacement weights not emitted'
old_volume=V(old['nodes'][old['tets']]);parents=ConstrainedCage(old['nodes'],old['tets'],old['classes']).embed(cage_tets.mean(1))[2]
assert np.max(np.abs(np.bincount(parents,weights=cell_volume,minlength=len(old_volume))-old_volume))<1e-8*old_volume.max()
regularizer=original['materialQuadrature']['emptyElementRegularizer'];normalization=1+regularizer*len(old_volume)
background=regularizer/normalization*cell_volume/old_volume[parents]
weights=physical/total/normalization+background
assert np.all(weights>0)and abs(weights.sum()-1)<1e-8
# Check that preservation of background alone reconstructs the old stored
# refinement density formula; geometry/parent containment stay independent.
assert np.max(np.abs(a['weights']-old['weights'][parents]*cell_volume/old_volume[parents]))<1e-10
np.savez_compressed(out/'cage.npz',**{**{k:a[k]for k in a.files},'weights':weights,'sourceMaterialVolumes':physical,'backgroundWeights':background})
for name in ['conforming.npz','conforming.json','targets.json']:
 if(root/name).exists():shutil.copy2(root/name,out/name)
newmeta={**meta,'materialQuadrature':{'method':report['method'],'sourceVolumeMm3':total,'physicalCells':int((physical>0).sum()),'backgroundDensityPreserved':True,'backgroundWeightTotal':float(background.sum()),'originalRegularizerPerParent':regularizer},'previousIntegrationInputs':report['inputsSha256']}
(out/'preparation.json').write_text(json.dumps(newmeta,indent=2)+'\n')
report['weightsEmitted']=True;report['backgroundDensityPreserved']=True;report['newWeightsSha256']=hashlib.sha256((out/'cage.npz').read_bytes()).hexdigest();(out/'integration.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps(report),flush=True)
