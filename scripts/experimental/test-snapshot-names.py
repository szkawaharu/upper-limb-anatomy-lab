from pathlib import Path
from snapshot_names import snapshot_name
cases=[(3,0,30),(3.375,0,30),(3.75,0,30),(3.375,90,30),(3.375,0,3000),(3.3750000000001,0,30),(-3.375,0,30)]
names=[snapshot_name(*x)for x in cases];assert len(set(names))==len(cases)
for name in names:assert Path(name).with_suffix('.json').stem==name
assert names[0]=='reference-e3-p0-k30'
try:snapshot_name(float('nan'),0,30)
except ValueError:pass
else:raise AssertionError('Nonfinite snapshot name accepted')
print('PASS snapshot persistence: fractional elevations/planes/material values have distinct suffix-safe keys; no silent truncation or overwrite')
