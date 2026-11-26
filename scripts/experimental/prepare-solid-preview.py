"""Export one completed offline reference and the identical source binding.
No geometry decimation; the browser test never loads this into production.
"""
import sys,json,gzip,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
binding_file=Path(sys.argv[1]);npz=Path(sys.argv[2]);report_file=Path(sys.argv[3]);out=Path(sys.argv[4]);out.parent.mkdir(parents=True,exist_ok=True)
report=json.loads(report_file.read_text())
if not report.get('converged'):raise ValueError('Preview requires a converged full-space reference')
B=next(b for b in json.loads(binding_file.read_text())if b['spec']['id']=='v60DeltMid');a=np.load(npz);n=len(B['rest'])//3
assert np.max(np.abs(a['referenceWorld'][:n].ravel()-B['rest']))<1e-12,'Source rest points differ'
assert np.array_equal(a['sourceFaces'].ravel(),B['faces']),'Source faces differ'
assert np.array_equal(a['fixed'],np.flatnonzero(np.isin(B['weights'],[0,1]))),'Source footprint masks differ'
assert np.all(np.isfinite(a['world']))
targets_file=npz.parent/'targets.json';target=next(x for x in json.loads(targets_file.read_text())['poses']if x['elev']==0 and x['plane']==0)
rest=np.array(B['rest']).reshape(-1,3);weights=np.array(B['weights'])
for key,weight in [('a',0),('b',1)]:
 bone=target[key];ids=np.flatnonzero(weights==weight);expected=Rotation.from_quat(bone['q']).apply(rest[ids]-bone['rest'])+bone['p']
 assert np.max(np.abs(expected-a['world'][ids]))<1e-10,'Saved solution does not have the declared neutral boundary pose'
packet={'experimental':True,'accepted':False,'id':'v60DeltMid','pose':{'elev':0,'plane':0},'target':target,'binding':B,'positions':a['world'][:n].ravel().tolist(),'report':report,'inputs':{str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in [binding_file,npz,report_file,targets_file]}}
with gzip.open(out,'wt',encoding='utf-8')as f:json.dump(packet,f,separators=(',',':'))
print(json.dumps({'file':str(out),'sourceVertices':n,'sourceFaces':len(B['faces'])//3,'experimental':True,'accepted':False}))
