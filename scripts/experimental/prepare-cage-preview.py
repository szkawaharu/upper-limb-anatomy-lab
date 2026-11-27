"""Export a fully converged conforming-cage snapshot for test-only rendering.
Rejects partial load paths. The source surface and its refinement have distinct
face arrays; their rest geometry is identical, their triangulations are not.
"""
import sys,json,gzip,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation

binding_file,npz,out=map(Path,sys.argv[1:4]);report_file=npz.with_suffix('.json')
report=json.loads(report_file.read_text());a=np.load(npz)
if 'snapshotSha256'in report:assert hashlib.sha256(npz.read_bytes()).hexdigest()==report['snapshotSha256'],'Mixed or incomplete reference snapshot'
if 'progress'in a:assert float(a['progress'])==report['progress']
assert report['progress']==1 and report['allStagesConverged']
assert all(s['converged']for s in report.get('stages', []))
assert all(s['minimization']['converged']for s in report.get('history', []))
B=next(b for b in json.loads(binding_file.read_text())if b['spec']['id']=='v60DeltMid')
candidate_file=Path(sys.argv[4])if len(sys.argv)>4 else binding_file
candidate=next(b for b in json.loads(candidate_file.read_text())if b['spec']['id']=='v60DeltMid')
rest=np.asarray(B['rest']).reshape(-1,3);n=len(rest)
assert np.max(np.abs(a['referenceWorld'][:n]-rest))<1e-12
assert np.array_equal(a['sourceWeights'][:n],candidate['weights'])
assert np.array_equal(B['rest'],candidate['rest'])and np.array_equal(B['faces'],candidate['faces'])
assert np.isfinite(a['sourceWorld']).all()
for key,weight in [('a',0),('b',1)]:
 bone=report['target'][key];ids=np.flatnonzero(a['sourceWeights']==weight)
 expected=Rotation.from_quat(bone['q']).apply(a['referenceWorld'][ids]-bone['rest'])+bone['p']
 assert np.max(np.abs(expected-a['sourceWorld'][ids]))<1e-10
packet={'experimental':True,'accepted':False,'kind':'conforming-cage','id':'v60DeltMid',
 'pose':{key:report['target'][key]for key in ['elev','plane']},'target':report['target'],'binding':B,
 'surface':{'rest':a['referenceWorld'].ravel().tolist(),'faces':a['sourceFaces'].ravel().tolist(),'weights':a['sourceWeights'].tolist()},
 'positions':a['sourceWorld'].ravel().tolist(),'report':{'converged':True,'volumeRatio':report['surfaceVolumeRatio'],'reference':report},
 'boundaryComparison':{'identical':bool(np.array_equal(B['weights'],candidate['weights'])),'baselineOriginVertices':int((np.array(B['weights'])==0).sum()),'candidateOriginVertices':int((np.array(candidate['weights'])==0).sum()),'candidateIsMeasuredEnthesis':False},
 'inputs':{str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in [binding_file,candidate_file,npz,report_file]}}
out.parent.mkdir(parents=True,exist_ok=True)
with gzip.open(out,'wt',encoding='utf-8')as f:json.dump(packet,f,separators=(',',':'))
print(json.dumps({'file':str(out),'sourceVertices':n,'displayVertices':len(a['referenceWorld']),'pose':packet['pose'],'accepted':False}))
