"""Retain compact, source-hashed evidence; never turn a partial run into pass."""
import sys,json,hashlib
from pathlib import Path
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
for arg in sys.argv[2:]:
 source=Path(arg);raw=json.loads(source.read_text());data={k:v for k,v in raw.items()if k not in ['rows','stages']}
 if 'rows'in raw:data['iterationOrStageCount']=len(raw['rows']);data['lastRow']=raw['rows'][-1]if raw['rows']else None
 if 'stages'in raw:data['stages']=[{k:v for k,v in r.items()if k!='rows'}for r in raw['stages']]
 evidence={'experimental':True,'accepted':False,'evidenceSource':str(source),'evidenceSha256':hashlib.sha256(source.read_bytes()).hexdigest(),'data':data}
 name=source.parent.name+'--'+source.name;(out/name).write_text(json.dumps(evidence,indent=2)+'\n');print(name)
