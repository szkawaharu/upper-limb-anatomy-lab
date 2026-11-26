"""Read pinned STL coordinates; optional offline atlas audit (numpy/scipy/matplotlib)."""
import struct, json, sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import least_squares
import hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(sys.argv[1] if len(sys.argv)>1 else '/tmp/upper-limb-atlas')
out=Path(sys.argv[2] if len(sys.argv)>2 else 'test-results/atlas-audit');out.mkdir(parents=True,exist_ok=True)
def stl(name):
    b=(root/name).read_bytes();n=struct.unpack_from('<I',b,80)[0]
    dtype=np.dtype([('n','<f4',(3,)),('v','<f4',(3,3)),('a','<u2')])
    v=np.frombuffer(b,offset=84,count=n,dtype=dtype)['v'].reshape(-1,3).astype(float)
    return np.unique(v,axis=0)
bones={k:stl('FMA'+str(i)+'.stl') for k,i in {'clavicle':13322,'scapula':13395,'humerus':23130,'manubrium':7486}.items()}
old={'SC':[-15,-80,1343],'AC':[-137,-92,1342],'GH':[-172,-72,1320]}
fig,axs=plt.subplots(1,3,figsize=(16,7))
for ax,(i,j,name) in zip(axs,[(0,2,'coronal x-z'),(0,1,'transverse x-y'),(1,2,'sagittal y-z')]):
    for k,p in bones.items():ax.scatter(p[::max(1,len(p)//6000),i],p[::max(1,len(p)//6000),j],s=.3,label=k)
    for k,p in old.items():ax.plot(p[i],p[j],'rx');ax.annotate('old '+k,(p[i],p[j]))
    ax.set_aspect('equal');ax.set_title(name+' (atlas mm)');ax.grid(alpha=.2)
axs[0].legend();fig.tight_layout();fig.savefig(out/'bone-landmark-audit.png',dpi=140)
c=bones['clavicle'];s=bones['scapula'];m=bones['manubrium']
# The articular candidate patches are terminal parts of this mesh, not the
# mesh origin; 8mm selects the terminal caps and is audited below.
ends={'SC':c[c[:,0]>c[:,0].max()-8],'AC':c[c[:,0]<c[:,0].min()+8]}
result={'old':old,'bounds':{k:[p.min(0).tolist(),p.max(0).tolist()] for k,p in bones.items()},'pairs':{}}
for name,patch in ends.items():
    other=m if name=='SC' else s
    dist,index=cKDTree(other).query(patch)
    ix=dist<=dist.min()+1.5
    pair=(patch[ix]+other[index[ix]])/2
    result['pairs'][name]={'clavicularCapCenter':patch.mean(0).tolist(),'closestPatchCenter':pair.mean(0).tolist(),'clavicularPatchCenter':patch[ix].mean(0).tolist(),'counterpartPatchCenter':other[index[ix]].mean(0).tolist(),'minGapMm':float(dist.min()),'maxPatchGapMm':float(dist[ix].max()),'n':int(ix.sum())}
# Candidate head-only sphere: medial-facing proximal articular cap.
h=bones['humerus'];cap=h[(h[:,2]>h[:,2].max()-40)&(h[:,0]>-172)]
A=np.c_[2*cap,np.ones(len(cap))];fit=np.linalg.lstsq(A,np.sum(cap**2,axis=1),rcond=None)[0];center=fit[:3];radius=np.sqrt(fit[3]+center@center)
res=np.linalg.norm(cap-center,axis=1)-radius
result['humeralSphere']={'center':center.tolist(),'radius':float(radius),'rms':float(np.sqrt(np.mean(res**2))),'max':float(abs(res).max()),'n':len(cap),'method':'proximal-medial cap least-squares; manual region estimate; not validated GH'}
(out/'landmarks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))

# Fixed, documented source regions identify distal condylar circular sections.
# This estimates a geometric hinge; a static atlas cannot measure a functional IHA.
ends=[]; fits=[]
for x in [-230,-214]:
    p=h[(abs(h[:,0]-x)<3)&(h[:,2]<1065)]
    f=least_squares(lambda a:np.linalg.norm(p[:,1:3]-a[:2],axis=1)-a[2],[-73,1058,12],loss='soft_l1',f_scale=1)
    res=np.linalg.norm(p[:,1:3]-f.x[:2],axis=1)-f.x[2]
    ends.append([float(p[:,0].mean()),*map(float,f.x[:2])]);fits.append({'radiusMm':float(f.x[2]),'rmsMm':float(np.sqrt(np.mean(res**2))),'n':len(p)})
result['elbow']={'axisEnds':ends,'center':np.mean(ends,axis=0).tolist(),'fits':fits,'method':'source distal condylar circular sections, x-band3mm, z<1065mm, robust soft-l1'}
cap=h[(h[:,2]>1305)&(h[:,0]>-172)]
f=least_squares(lambda a:np.linalg.norm(cap-a[:3],axis=1)-a[3],[-162.69,-72.24,1325.72,21.5],loss='soft_l1',f_scale=1)
result['robustHumeralSphere']={'center':f.x[:3].tolist(),'radiusMm':float(f.x[3]),'rmsMm':float(np.sqrt(np.mean((np.linalg.norm(cap-f.x[:3],axis=1)-f.x[3])**2))),'n':len(cap),'method':'proximal medial region z>1305,x>-172; robust soft-l1'}
result['sources']={file:hashlib.sha256((root/file).read_bytes()).hexdigest() for file in ['FMA13322.stl','FMA13395.stl','FMA23130.stl','FMA7486.stl']}
result['atlasCommit']='2c83c698788cff57f6317bddcf712e7acfbea7d8'
(out/'landmarks.json').write_text(json.dumps(result,indent=2)+'\n')

# ISB thorax landmark estimates from the same pinned atlas.
landmarks={}
for name,file in [('C7','FMA12525.stl'),('T8','FMA9991.stl'),('PX','FMA7488.stl'),('IJ','FMA7486.stl')]:
    p=stl(file)
    if name in ['C7','T8']: s=p[p[:,1]>p[:,1].max()-3]
    elif name=='PX': s=p[p[:,2]<p[:,2].min()+2]
    else:
        s=p[abs(p[:,0])<5];s=s[s[:,2]>s[:,2].max()-2]
    landmarks[name]={'point':s.mean(0).tolist(),'samples':len(s),'file':file,'sha256':hashlib.sha256((root/file).read_bytes()).hexdigest(),'confidence':'source-surface extremum estimate, not digitized palpation'}
result['thoraxLandmarks']=landmarks
(out/'landmarks.json').write_text(json.dumps(result,indent=2)+'\n')

# Overlay the new source estimates on the unchanged bone cloud. Labels are
# deliberately geometric estimates, not clinically measured joint centres.
new={'SC':result['pairs']['SC']['clavicularCapCenter'],'AC':result['pairs']['AC']['closestPatchCenter'],'GH':result['robustHumeralSphere']['center']}
for ax,(i,j,_) in zip(axs,[(0,2,''),(0,1,''),(1,2,'')]):
    for name,p in new.items():
        ax.plot(p[i],p[j],'g+',markersize=10)
        ax.annotate('estimated '+name,(p[i],p[j]),xytext=(4,-12),textcoords='offset points',fontsize=8)
    ax.set_xlabel(['x','y','z'][i]+' / mm');ax.set_ylabel(['x','y','z'][j]+' / mm')
fig.tight_layout();fig.savefig(out/'bone-landmark-audit.png',dpi=140)
