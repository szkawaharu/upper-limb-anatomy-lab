"""Plot measured saved-path audits, including the rejected local-J outcome."""
import sys,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root,out=map(Path,sys.argv[1:3]);path=json.loads((root/'stopped-path.json').read_text());audit=json.loads((root/'independent-audit.json').read_text());rows=audit['rows'];angle=[r['elev']for r in rows]
fig,ax=plt.subplots(2,2,figsize=(11,8))
ax[0,0].plot(angle,[r['sourceVolumeRatio']for r in rows]);ax[0,0].axhline(1,c='gray',ls='--');ax[0,0].set_ylabel('Source surface volume ratio');ax[0,0].ticklabel_format(useOffset=False)
for key,index,label in [('surfaceAdjacentJRange',0,'Source-adjacent minimum'),('surfaceAdjacentJRange',1,'Source-adjacent maximum'),('cageJRange',1,'All cage maximum (incl. background)')]:ax[0,1].semilogy(angle,[r[key][index]for r in rows],label=label)
ax[0,1].set_ylabel('Local cell volume ratio J (log scale)');ax[0,1].axhline(1,c='gray',ls='--');ax[0,1].legend(fontsize=8)
ax[1,0].plot(angle,[r['maxDepthIncreaseMm']for r in rows]);ax[1,0].axhline(.1,c='gray',ls='--',label='Sampled contact stopping target');ax[1,0].set_ylabel('Maximum new bone depth at vertices (mm)');ax[1,0].legend(fontsize=8)
sequence=path['rows'];previous=90.;degrees=[];movement=[]
for r in sequence:
 degrees.append(abs(r['elev']-previous));movement.append(r['maxStepWorldMm']);previous=r['elev']
scatter=ax[1,1].scatter(angle,movement,c=degrees,cmap='viridis');fig.colorbar(scatter,ax=ax[1,1],label='Actual accepted angle step (degrees)');ax[1,1].set_ylabel('Maximum source-vertex step in world (mm)')
for a in ax.flat:a.set_xlabel('Humerothoracic abduction (degrees)');a.grid(alpha=.25)
fig.suptitle('Actual-angle contact reference / stopped at 147.375 degrees\nGlobal volume and contact do not qualify local deformation',fontsize=13)
fig.text(.5,.015,'Right middle deltoid only, neck zero. Old integration weights; local J remains unacceptable.\n55 converged poses; 180 degrees and the reverse path were not reached. Not deployed, no mobile/GPU solver certification.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.07,1,.93]);out.parent.mkdir(parents=True,exist_ok=True);fig.savefig(out,dpi=150);svg=out.with_suffix('.svg');fig.savefig(svg);plt.close(fig);svg.write_text('\n'.join(line.rstrip()for line in svg.read_text().splitlines())+'\n');print(out)
