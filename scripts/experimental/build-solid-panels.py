"""Scientific comparison panels from unchanged browser screenshots.
No geometry/image synthesis or retouching. Raw images remain beside panels.
"""
import sys,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
p=Path(sys.argv[1]);data=json.loads((p/'comparison.json').read_text());assert len(data['rows'])==12 and data['errors']==[]
for isolated in [False,True]:
 fig,axes=plt.subplots(3,2,figsize=(10,13),facecolor='white')
 for row,view in enumerate(['front','side','back']):
  pair=[next(x for x in data['rows']if x['view']==view and x['mode']==mode and x['isolation']==isolated)for mode in ['curve','continuum']]
  assert pair[0]['camera']==pair[1]['camera']
  for col,item in enumerate(pair):
   axes[row,col].imshow(mpimg.imread(p/item['file']));axes[row,col].axis('off');axes[row,col].set_title(item['mode']+' / '+view,fontsize=11)
 fig.suptitle(str(data['pose'])+' / '+('isolated deltoid middle'if isolated else'all muscles')+'\nIdentical rest surface and final bones; experimental, not accepted',fontsize=12)
 fig.subplots_adjust(top=.93,bottom=.01,left=.01,right=.99,wspace=.02,hspace=.08)
 file=p/('panel-'+('isolated'if isolated else'all')+'.png');fig.savefig(file,dpi=120);plt.close(fig);print(file)
