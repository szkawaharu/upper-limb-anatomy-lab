"""Diagnostic only: source fixed attachment occupancy vs GH axial rotation.
Scapular axes use TS->AC. The exploratory Y-X-Y rotation sequence has NOT
been matched to the Ludewig2009 reporting sequence. Trial angles are sensitivity
inputs, not adopted physiological targets. No production rig is mutated.
"""
import json
from pathlib import Path
import numpy as np,trimesh
from scipy.spatial.transform import Rotation
p=Path('test-results/v076-cage-mid-18-projected');a=np.load(p/'cage.npz');raw=np.load('test-results/v076-volumetric-quality/deltoid-mid-source.npz');cal=json.loads(Path('docs/validation/v076/source-bones/browser-numeric/rig-calibration.json').read_text());L={k:np.array(v)for k,v in cal['landmarks'].items()};scale=cal['scale'];P=np.array([[-1,0,0],[0,0,1],[0,1,0]],float);offset=L['SC']/scale-P@np.array([-14.91024,-129.66467,1341.37453]);mesh=trimesh.load_mesh('/tmp/upper-limb-atlas/FMA23130.stl',process=True);mesh.vertices=np.asarray(mesh.vertices)@P.T+offset
norm=lambda x:x/np.linalg.norm(x)
z=norm(L['AC']-L['TS']);x=norm(np.cross(L['TS']-L['AC'],L['AI']-L['AC']))
if x@[0,0,-1]<0:x=-x
S0=np.column_stack([x,norm(np.cross(z,x)),z]);y=norm(L['GH']-L['EL']);x=norm(np.cross(y,L['EL_LAT']-L['EL_MED']));z=norm(np.cross(x,y));H0=np.column_stack([x,y,z]);rows=[]
for plane in [0,90]:
 target=json.loads((p/('reference-e90-p'+str(plane)+'-k30.json')).read_text())['target'];A=target['a'];B=target['b'];RA=Rotation.from_quat(A['q']).as_matrix();RB=Rotation.from_quat(B['q']).as_matrix();S=RA@S0;relative=S.T@RB@H0;v=relative[:,1];beta=-np.arccos(np.clip(v[1],-1,1));alpha=np.arctan2(-v[0],-v[2]);noax=Rotation.from_rotvec([0,alpha,0]).as_matrix()@Rotation.from_rotvec([beta,0,0]).as_matrix();ax=noax.T@relative;current=np.rad2deg(np.arctan2(ax[0,2],ax[0,0]));points=raw['p'][raw['weights']==0];posed=(points-A['rest'])@RA.T+A['p'];row={'plane':plane,'elev':90,'GHCurrentAxialDeg':float(current),'GHElevationDeg':float(np.rad2deg(beta)),'GHPlaneDeg':float(np.rad2deg(alpha)),'fixedOriginVertices':len(points),'tests':[]}
 for gamma in [current,-45,-60,-75]:
  rotated=S@noax@Rotation.from_rotvec([0,np.deg2rad(gamma),0]).as_matrix()@H0.T;assert np.linalg.norm(rotated@y-RB@y)<1e-12
  query=((posed-B['p'])@rotated+np.array(B['rest']))/scale;inside=mesh.contains(query);depth=trimesh.proximity.closest_point(mesh,query[inside])[1]if inside.any()else np.array([])
  row['tests'].append({'GHReportAxialDeg':float(gamma),'originInsideHumerus':int(inside.sum()),'maxDepthMm':float(depth.max(initial=0)),'humerusQ':Rotation.from_matrix(rotated).as_quat().tolist(),'humeralShaftChange':float(np.linalg.norm(rotated@y-RB@y))})
 rows.append(row);print(json.dumps(row),flush=True)
(p/'axial-diagnostic.json').write_text(json.dumps({'experimental':True,'accepted':False,'method':__doc__,'rows':rows},indent=2)+'\n')
