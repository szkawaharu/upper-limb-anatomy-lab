"""Independent per-source-element coverage audit; never qualifies anatomy."""
import sys,json,hashlib
from pathlib import Path
import numpy as np
root,fine_path,out=map(Path,sys.argv[1:4]);a=np.load(root/'cage.npz');fine=np.load(fine_path);checkpoint=root/'integration-progress.npz';coverage=np.load(checkpoint);assert int(coverage['cellsCompleted'])==len(a['tets'])
points=(fine['nodes']-a['centerWorld'])/float(a['worldUnitsPerMm']);p=points[fine['tets']]
# Independent scalar triple product, not the integration driver's determinant.
volume=np.einsum('ni,ni->n',np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]),p[:,3]-p[:,0])/6
assert np.all(volume>0)
covered=coverage['coveredFineVolumes'];error=np.abs(covered-volume);relative=error/volume;idx=int(relative.argmax());zero=covered==0
r={'experimental':True,'accepted':False,'sourceFineTets':len(volume),'relativeErrorQuantiles':np.quantile(relative,[0,.5,.99,.999,1]).tolist(),'worstRelativeFineTet':{'id':idx,'relativeError':float(relative[idx]),'sourceVolumeMm3':float(volume[idx]),'absoluteErrorMm3':float(error[idx])},'fineTetsWithZeroReportedIntersection':int(zero.sum()),'maximumZeroCoverageFineTetVolumeMm3':float(volume[zero].max(initial=0)),'sumAbsoluteCoverageErrorMm3':float(error.sum()),'inputSha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in [root/'cage.npz',fine_path,checkpoint]},'codeSha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'limitations':['Floating clipping may be uncertain for near-contact slivers; not exact rational geometry','Small global error alone is not local or anatomical acceptance']}
out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
