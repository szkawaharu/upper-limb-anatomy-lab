"""One source-conforming, stateless continuum reference pose. Offline only.
Boundary continuation is a numerical solve path, not measured muscle activity.
"""
import json,time,sys,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from volumetric_reference import Solid

directory=Path(sys.argv[1]if len(sys.argv)>1 else'test-results/v076-volumetric')
elev=float(sys.argv[2]if len(sys.argv)>2 else 0);plane=float(sys.argv[3]if len(sys.argv)>3 else 0)
iterations=int(sys.argv[4]if len(sys.argv)>4 else 200)
initial_bulk=float(sys.argv[5]if len(sys.argv)>5 else 30)
run_tag=sys.argv[6]if len(sys.argv)>6 else'v2'
targetData=json.loads((directory/'targets.json').read_text());target=next(x for x in targetData['poses']if x['elev']==elev and x['plane']==plane)
if targetData['id']!='v60DeltMid':raise ValueError('This first reference driver requires the audited deltoid-middle fixture')
source=np.load(directory/'deltoid-mid-source.npz');mesh=np.load(directory/'deltoid-mid-tetrahedra.npz');world=mesh['nodes'];cells=mesh['tets'];weights=source['weights'];scale=targetData['worldUnitsPerMm'];center=world.mean(0);rest=(world-center)/scale
fixed=np.flatnonzero((weights==0)|(weights==1));classes=weights[fixed].astype(int);start=time.time();solid=Solid(rest,cells,fixed)
assert np.array_equal(source['p'],world[:len(weights)]),'Source boundary vertex order changed'
assert np.all(solid.inactive>=len(weights)),'An original anatomical vertex was lost from the tetrahedral mesh'
def axis_angle(bone):
 q=np.array(bone['q']);s=np.linalg.norm(q[:3]);angle=2*np.arctan2(s,q[3])
 if s<1e-12:
  if q[3]<0:raise ValueError('Ambiguous full-turn quaternion lift')
  return np.zeros(3),0.
 return q[:3]/s,float(angle)
angles=[axis_angle(target[k])for k in ['a','b']];max_angle=max(a[1]for a in angles)
def pins(progress):
 result=np.zeros((len(fixed),3))
 for k,key in enumerate(['a','b']):
  ids=np.flatnonzero(classes==k);bone=target[key];axis,angle=angles[k];q=Rotation.from_rotvec(axis*angle*progress);anchor=np.array(bone['rest']);at=anchor+progress*(np.array(bone['p'])-anchor);points=q.apply(world[fixed[ids]]-anchor)+at;result[ids]=(points-center)/scale
 return result
x=rest.copy();progress=0.;base_step=min(1.,np.deg2rad(3)/max(max_angle,1e-12));step=base_step;rows=[];bulk=initial_bulk
while progress<1-1e-12:
 next_progress=min(1.,progress+step);predicted=solid.predict_boundary(x,pins(next_progress))
 if not np.isfinite(solid.evaluate(predicted,bulk,gradient=False)[0]):
  step*=.5
  if step<1e-7:raise RuntimeError('Could not advance positive-J Dirichlet continuation')
  continue
 x,report=solid.minimize(predicted,bulk=bulk,iterations=iterations)
 while abs(report['volumeRatio']-1)>.05 and bulk<3000:
  bulk*=10;x,report=solid.minimize(x,bulk=bulk,iterations=iterations)
 progress=next_progress;row={'progress':progress,'step':step,'elapsedSeconds':time.time()-start,**report};rows.append(row);print(json.dumps(row),flush=True)
 path=directory/('reference-'+run_tag+'-e'+str(int(elev))+'-p'+str(int(plane))+'-i'+str(iterations)+'-k'+str(int(initial_bulk)));np.savez(str(path)+'.npz',world=x*scale+center,tets=cells,referenceWorld=world,sourceFaces=source['faces'],fixed=fixed)
 result={'experimental':True,'productionEnabled':False,'model':'normalized isochoric distortion plus log-volume penalty; no measured material constants','input':targetData,'target':target,'iterationsPerStage':iterations,'rows':rows,'completeBoundaryProgress':progress,'allStagesConverged':all(r['converged']for r in rows),'attachmentMaxMm':float(np.linalg.norm(x[fixed]-pins(progress),axis=1).max()),'sourceBoundaryPreserved':True,'unreferencedGeneratedNodes':solid.inactive.tolist(),'boneContactEvaluated':False,'elapsedSeconds':time.time()-start,'codeSha256':{p:hashlib.sha256(Path(__file__).with_name(p).read_bytes()).hexdigest()for p in ['volumetric_reference.py','run-volumetric-reference.py']}}
 Path(str(path)+'.json').write_text(json.dumps(result,indent=2)+'\n');step=min(base_step,step*1.25)
print('FINISHED',json.dumps({'elev':elev,'plane':plane,'stages':len(rows),'seconds':time.time()-start,'last':rows[-1]}),flush=True)
