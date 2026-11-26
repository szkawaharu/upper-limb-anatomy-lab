"""Exact source-face volume integral. This is signed enclosed volume, not a
local incompressibility, non-inversion or collision certificate. No thresholds
are adjusted. Rigid fixed faces are integrated analytically, every mixed
source triangle is retained (multiplier one).
"""
import json,sys
from pathlib import Path
import numpy as np
src=Path(sys.argv[1]if len(sys.argv)>1 else'test-results/v076-source-graphs/bindings.json');out=Path(sys.argv[2]if len(sys.argv)>2 else'test-results/v076-source-graphs/volume');out.mkdir(parents=True,exist_ok=True)
for B in json.loads(src.read_text()):
 if not B['meshPreparation']['closed']or not B['meshPreparation']['oriented']:continue
 p=np.array(B['rest']).reshape(-1,3);f=np.array(B['faces']).reshape(-1,3);w=np.array(B['weights']);origin=p.mean(0);P=p[f]-origin;control=np.einsum('ij,ij->i',P[:,0],np.cross(P[:,1],P[:,2]));reference=float(control.sum()/6)
 varying=B['spec']['id']in['v58Levator','v59TrapUpper'];pureA=(w[f]==0).all(1)if not varying else np.zeros(len(f),bool);pureB=(w[f]==1).all(1);mixed=~(pureA|pureB);F=f[mixed];used=np.unique(F);remap=np.full(len(p),-1,int);remap[used]=np.arange(len(used));fn=np.cross(p[f[:,1]]-p[f[:,0]],p[f[:,2]]-p[f[:,0]])
 Q={'sourceVertices':used.tolist(),'faces':remap[F].ravel().tolist(),'multipliers':np.ones(len(F)).tolist(),'controlDet':control[mixed].tolist(),'referenceVolume':reference,'sourceCenter':origin.tolist(),'pureBVector':fn[pureB].sum(0).tolist(),'method':__doc__,'sourceMixedFaces':len(F)}
 (out/(B['spec']['id']+'-quadrature.json')).write_text(json.dumps(Q,separators=(',',':')));print(B['spec']['id'],len(F),len(used),flush=True)
