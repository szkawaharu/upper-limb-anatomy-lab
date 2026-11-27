"""Compare original/corrected material integration at identical scene settings."""
import sys,json,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
old,new,out=map(Path,sys.argv[1:4]);out.mkdir(parents=True,exist_ok=True)
a=json.loads((old/'comparison.json').read_text());b=json.loads((new/'comparison.json').read_text());assert a['pose']==b['pose']
for isolation in [False,True]:
 fig,axes=plt.subplots(3,2,figsize=(10,13))
 for row,view in enumerate(['front','side','back']):
  left=next(r for r in a['rows']if r['view']==view and r['isolation']==isolation and r['mode']=='continuum');right=next(r for r in b['rows']if r['view']==view and r['isolation']==isolation and r['mode']=='continuum');assert left['camera']==right['camera']
  for bone in ['a','b']:
   for key in ['q','p','rest']:assert np.max(np.abs(np.array(left['audit'][bone][key])-right['audit'][bone][key]))<1e-12
  for col,(root,data,label)in enumerate([(old,left,'Original integration'),(new,right,'Source-intersection integration')]):
   axes[row,col].imshow(plt.imread(root/data['file']));axes[row,col].axis('off');axes[row,col].set_title(label+' / '+view,fontsize=11)
 fig.suptitle('90-degree abduction / '+('isolated middle deltoid'if isolation else'all muscles')+'\nSame rest mesh, attachments, final bones and cameras; local J still fails',fontsize=12)
 fig.tight_layout(rect=[0,0,1,.95]);path=out/('panel-integration-'+('isolated'if isolation else'all')+'.png');fig.savefig(path,dpi=130);plt.close(fig);print(path)
(out/'manifest.json').write_text(json.dumps({'experimental':True,'accepted':False,'sameFinalBonesAndCameras':True,'inputsSha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in [old/'comparison.json',new/'comparison.json']},'limitations':['One reference muscle/pose','Other muscles remain production curve engine','No medical, whole-ROM or mobile/GPU solver acceptance']},indent=2)+'\n')
