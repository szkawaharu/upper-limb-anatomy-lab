import json,re
from pathlib import Path
from bone_sets import bone_set,DEFINITIONS
index=Path('index.html').read_text();data=json.loads(DEFINITIONS.read_text())
for id,entry in data['models'].items():
 match=re.search(r"\{id:'"+id+r"',[^\n]*a:'([^']+)',b:'([^']+)'",index);assert match,'Missing existing muscle definition '+id;roles=dict(zip(['a','b'],match.groups()))
 for bone in entry['bones']:assert roles[bone['key']]==bone['boneId'],'Reference bone role differs from application attachment'
 assert entry['limitations'],'Partial contact coverage must be explicit'
try:bone_set('unregistered-muscle')
except ValueError:pass
else:raise AssertionError('Unknown muscle silently used another bone geometry')
print('PASS bone adapter integration: body roles match existing application definitions, partial coverage explicit, unknown models rejected')
