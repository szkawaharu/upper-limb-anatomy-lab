"""Production clustered trace mesh vs independently repaired full source.
The different rest preparations and estimated boundary data are explicit.
"""
import sys,json,gzip,hashlib
from pathlib import Path
import numpy as np
binding,npz,tracefile,out=map(Path,sys.argv[1:5]);trace=json.loads(tracefile.read_text());B=next(b for b in json.loads(binding.read_text())if b['spec']['id']==trace['id']);a=np.load(npz);reportfile=npz.with_suffix('.json');r=json.loads(reportfile.read_text());assert r['progress']==1 and r['allStagesConverged'];assert hashlib.sha256(npz.read_bytes()).hexdigest()==r['snapshotSha256'];expected=next(row for row in trace['rows']if row['elev']==r['target']['elev']);packet={'experimental':True,'accepted':False,'kind':'full-source-reference','id':trace['id'],'restSurfaceIdentical':False,'pose':{'elev':r['target']['elev'],'plane':r['target']['plane']},'target':r['target'],'binding':B,'surface':{'rest':a['referenceWorld'].ravel().tolist(),'faces':a['sourceFaces'].ravel().tolist(),'weights':a['sourceWeights'].tolist()},'positions':a['sourceWorld'].ravel().tolist(),'traceExpected':{'vertex':trace['vertex'],'point':expected['point']},'report':{'converged':True,'volumeRatio':r['surfaceVolumeRatio'],'reference':r},'boundaryComparison':{'identical':False,'baseline':'captured production partial attachment field','candidate':'source-derived footprints including inferred insertion; still unverified'},'inputs':{str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in [binding,npz,reportfile,tracefile]}}
assert hashlib.sha256(binding.read_bytes()).hexdigest()==trace['inputSha256'],'Production trace binding changed'
out.parent.mkdir(parents=True,exist_ok=True)
with gzip.open(out,'wt')as stream:json.dump(packet,stream,separators=(',',':'))
print(json.dumps({'file':str(out),'id':trace['id'],'pose':packet['pose'],'restSurfaceIdentical':False,'productionTraceExpected':packet['traceExpected']}))
