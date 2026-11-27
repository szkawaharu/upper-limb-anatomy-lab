"""Explicit reference bone geometry selection; unknown models never fall back
to another muscle's scapula/humerus. Coordinate calculations remain in the
independent material-frame solver and world-frame auditor.
"""
import json
from pathlib import Path
DEFINITIONS=Path(__file__).with_name('contact-bone-sets.json')
def model_id(file,report):
 for key in ['contactModelId','modelId']:
  if key in report:return report[key]
 for p in [Path(file).parent,Path(file).parent.parent]:
  if(p/'targets.json').exists():return json.loads((p/'targets.json').read_text())['id']
 raise ValueError('Reference has no explicit model identity')
def bone_set(id):
 data=json.loads(DEFINITIONS.read_text())
 if id not in data['models']:raise ValueError('No audited bone geometry adapter for '+id)
 return data['models'][id]
