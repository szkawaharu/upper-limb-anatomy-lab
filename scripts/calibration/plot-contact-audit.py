"""Engineering kinematics and quadrature audit; no clinical classification."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path('docs/validation/v076/contact-rms')
data=json.loads((out/'shoulder-motion.json').read_text())
fig,axs=plt.subplots(1,3,figsize=(14,4.6),sharey=True)
for plane,ax in zip([0,40,90],axs):
    for model,label,color in [('before','before: sum','#777777'),('after','after: spatial mean','#007d79')]:
        rows=[r for r in data['rows'] if r['model']==model and r['input']['plane']==plane]
        ax.plot([r['input']['elev'] for r in rows],[r['motion']['scapulothoracic']['upwardRotation'] for r in rows],'o-',color=color,label=label)
    ax.plot([r['input']['elev'] for r in rows],[r['motion']['unconstrainedPriorScapulothoracic']['upwardRotation'] for r in rows],'--',color='#bd7b17',label='prior reconstruction')
    ax.axvspan(120,180,color='#eeeeee');ax.set_title(f'Plane {plane} deg');ax.set_xlabel('HT elevation input / degree');ax.grid(alpha=.2)
axs[0].set_ylabel('source-calibrated ST upward rotation / degree');axs[0].legend(fontsize=9)
fig.suptitle('Effect of contact sample aggregation on scapular kinematics')
fig.text(.5,.01,'Approximate model; the prior is reconstructed from SC/AC curves. Shaded high elevation is unvalidated. No disease diagnosis.',ha='center',fontsize=8)
fig.tight_layout(rect=[0,.04,1,.95]);fig.savefig(out/'scapular-motion.png',dpi=150);fig.savefig(out/'scapular-motion.svg')
fig,ax=plt.subplots(figsize=(8,4.5))
for model,label,color in [('before','before: sum','#777777'),('after','after: spatial mean','#007d79')]:
    rows=[r for r in data['density'] if r['model']==model and r['input']['plane']==0]
    ax.plot([1]+[r['copies'] for r in rows],[0]+[r['ACdifferenceMm'] for r in rows],'o-',label=label,color=color)
ax.set(xticks=[1,2,4],xlabel='copies of the SAME spatial samples',ylabel='AC position change / mm',title='Same geometry must not become a stronger constraint')
ax.grid(alpha=.2);ax.legend();fig.tight_layout();fig.savefig(out/'sampling-invariance.png',dpi=150)
