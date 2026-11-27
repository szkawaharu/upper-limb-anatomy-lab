"""Lossless finite scalar labels without a dot that Path.with_suffix strips."""
import math
def snapshot_name(elevation,plane,bulk):
 def label(value):
  value=float(value)
  if not math.isfinite(value):raise ValueError('Nonfinite snapshot coordinate')
  return (str(int(value))if value.is_integer()else repr(value)).replace('.','d').replace('-','m')
 return 'reference-e'+label(elevation)+'-p'+label(plane)+'-k'+label(bulk)
