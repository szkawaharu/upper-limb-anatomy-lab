"""Reproducible engineering plots; ratios are not anatomical normal values."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path('docs/validation/v076')
data=[json.loads((out/f).read_text()) for f in ['baseline-v075-report.json','accepted-report.json','rejected-projection-report.json']]
labels=['v0.75 baseline','v0.76 adopted','REJECTED projection']
x=np.arange(len(data[0]['rows']))
fig,axs=plt.subplots(1,2,figsize=(13,5))
for d,label,offset,color in zip(data,labels,[-.25,0,.25],['#777777','#007d79','#bd3b37']):
    for ax,metric in zip(axs,['edgeRatio','areaRatio']):
        ax.bar(x+offset,[max(a[metric][1] for a in r['audit']) for r in d['rows']],.24,label=label,color=color)
for ax,title in zip(axs,['Maximum edge length ratio','Maximum triangle area ratio']):
    ax.set_yscale('log');ax.set_title(title);ax.set_ylabel('posed / source rest (log scale)')
    ax.set_xticks(x,[f"{r['elev']} / {r['plane']}" for r in data[0]['rows']],rotation=35)
    ax.set_xlabel('UI elevation / plane, degrees');ax.grid(axis='y',alpha=.2);ax.axhline(1,color='black',lw=.7)
axs[0].legend(fontsize=9)
fig.suptitle('Priority six atlas regions: local shape quality remains UNACCEPTABLE',fontsize=12)
fig.text(.5,.01,'Same UI inputs; neutral calibration differs. Unclosed surfaces: volume not evaluable. No clinical accuracy claim.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.04,1,.95]);fig.savefig(out/'shape-outliers.png',dpi=150);fig.savefig(out/'shape-outliers.svg')
fig,ax=plt.subplots(figsize=(7,4.5))
for r in data[1]['refinement']:
    ax.plot([x['stepDeg'] for x in r['rows']],[x['maxMm'] for x in r['rows']],'o-',label=r['id']+' vertex '+str(r['vertexId']))
ax.set(xlabel='joint increment / degree',ylabel='maximum sampled displacement / mm',title='Step refinement: finite jumps removed at captured cases')
ax.legend();ax.grid(alpha=.25);fig.text(.5,.01,'Convergence does not establish natural muscle shape. Current upper trapezius remains excessive.',ha='center',fontsize=8)
fig.tight_layout(rect=[0,.04,1,1]);fig.savefig(out/'continuity-refinement.png',dpi=150)
