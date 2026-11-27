"""Plots actual frozen trace/reference data, including remaining J outliers.
No pose interpolation or image synthesis. Different mesh/boundary preparations
are stated on the figure; this is not a production acceptance certificate.
"""
import sys,json
from pathlib import Path
import numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
old=json.loads(Path(sys.argv[1]).read_text());new=json.loads(Path(sys.argv[2]).read_text());out=Path(sys.argv[3]);out.parent.mkdir(parents=True,exist_ok=True);scale=.00555;rows=new['rows'];angle=np.array([r['elev']for r in rows]);baseline=np.array([next(x['point']for x in old['rows']if x['elev']==v)for v in angle]);reference=np.array([r['point']for r in rows]);steps=new['steps'];mid=np.array([(s['from']+s['to'])/2 for s in steps]);fig,ax=plt.subplots(2,2,figsize=(11,8.5))
ax[0,0].plot(angle,np.linalg.norm(baseline-baseline[0],axis=1)/scale,'o-',label='Legacy curve / production trace 168');ax[0,0].plot(angle,np.linalg.norm(reference-reference[0],axis=1)/scale,'o-',label='Repaired-source reference / same spatial location');ax[0,0].set_ylabel('Distance from 3-degree position (mm)');ax[0,0].legend(fontsize=8)
ax[0,1].plot(mid,np.linalg.norm(np.diff(baseline,axis=0),axis=1)/scale,'o-',label='Legacy trace 168');ax[0,1].plot(mid,[s['sameSpatialLocationMm']for s in steps],'o-',label='Reference same location');ax[0,1].plot(mid,[s['sourceVertexMaxMm']for s in steps],'--',label='All reference vertices: maximum');ax[0,1].set_ylabel('Displacement per 0.375-degree step (mm)');ax[0,1].legend(fontsize=8)
ax[1,0].plot(angle,[r['volumeRatio']for r in rows],'o-');ax[1,0].axhline(1,c='gray',ls='--');ax[1,0].set_ylabel('Closed reference surface volume ratio');ax[1,0].ticklabel_format(useOffset=False)
ax[1,1].plot(angle,[r['JRange'][0]for r in rows],'o-',label='Minimum J');ax[1,1].plot(angle,[r['JRange'][-1]for r in rows],'o-',label='Maximum J');ax[1,1].axhline(1,c='gray',ls='--');ax[1,1].set_ylabel('Local element volume ratio J');ax[1,1].legend(fontsize=8)
for a in ax.flat:a.set_xlabel('Humerothoracic abduction (degrees)');a.grid(alpha=.25)
fig.suptitle('Upper trapezius / frozen final bones, neck zero\nDifferent mesh and attachment preparations; reference only, not deployed',fontsize=13)
fig.text(.5,.015,'Reference uses inferred attachment data and an independently repaired source mesh.\nLocal J outliers remain despite small global volume change. No anatomical, full-ROM, GPU-solver or mobile certification.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.07,1,.93]);fig.savefig(out,dpi=150);fig.savefig(out.with_suffix('.svg'));plt.close(fig);print(out)
