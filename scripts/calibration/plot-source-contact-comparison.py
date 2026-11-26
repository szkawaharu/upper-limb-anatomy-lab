"""Compose unchanged browser QA frames as labelled scientific comparison
panels. Every pair shares camera and UI pose inputs; no image content is edited.
"""
import json,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
src=Path(sys.argv[1]if len(sys.argv)>1 else'test-results/v076-source-bone-comparison');out=Path(sys.argv[2]if len(sys.argv)>2 else'test-results/v076-source-bone-comparison');out.mkdir(parents=True,exist_ok=True)
rows=json.loads((src/'comparison.json').read_text())
for name in dict.fromkeys(x['name']for x in rows):
 for mode in ['bones','muscles']:
  fig,axes=plt.subplots(3,2,figsize=(10,13.6),layout='constrained');fig.patch.set_facecolor('#ffffff')
  for i,view in enumerate(['front','side','back']):
   pair=[next(x for x in rows if x['name']==name and x['mode']==mode and x['view']==view and x['version']==version)for version in ['v075','v076']]
   assert pair[0]['camera']==pair[1]['camera']
   for j,version in enumerate(['v075','v076']):
    image=mpimg.imread(src/(version+'-'+name+'-'+mode+'-'+view+'.png'));ax=axes[i,j];ax.imshow(image);ax.set_axis_off();ax.set_title(version+' / '+view,fontsize=11)
  fig.suptitle(name+' / '+mode+'\nSame camera + UI angles; anatomy/volume quality remains unvalidated',fontsize=12)
  fig.savefig(out/('comparison-'+name+'-'+mode+'.png'),dpi=140);plt.close(fig)
