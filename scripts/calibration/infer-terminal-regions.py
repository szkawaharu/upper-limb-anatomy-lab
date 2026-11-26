"""Source-mesh terminal region candidates for absent Dirichlet ends.
Anatomical target bone sets are explicit. Only the prior end's half-domain is
searched; candidates are within one measured mesh edge spacing of its minimum
bone-surface distance. A muscle terminal with a tendon gap is NOT relabelled as
an observed enthesis. Existing exact attachment masks are retained unchanged.
"""
import json,re
from pathlib import Path
import numpy as np,trimesh
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
src=Path('test-results/v076-source-clean/bindings.json');out=Path('test-results/v076-source-terminals');out.mkdir(parents=True,exist_ok=True);cache=Path('/tmp/upper-limb-atlas')
fixture=json.load(open('docs/validation/v076/contact-rms/rig-calibration.json'));scale=fixture['scale'];SC=np.array([-14.91024,-129.66467,1341.37453]);P=np.array([[-1,0,0],[0,0,1],[0,1,0]],float);offset=np.array(fixture['landmarks']['SC'])/scale-P@SC
manifest=json.loads(re.search(r'const V57_MANIFEST=Object.freeze\((\[.*?\])\);',Path('index.html').read_text()).group(1));byid={x['id']:x['file']for x in manifest};meshes={}
def bone(key):
 if key not in meshes and key.startswith('scapula:'):
  base=bone('scapula');v=np.asarray(base.vertices);spacing=float(np.median(base.edges_unique_length));TS=np.array([-63.87020,-6.27689,1331.10999]);AA=np.array([-155.367,-51.18210,1356.91003]);AI=np.array([-91.16020,.02043,1201.5])
  if key.endswith('spine'):
   edges=base.edges_unique;normal=base.vertex_normals;start=int(np.argmin(np.linalg.norm(v-TS,axis=1)));end=int(np.argmin(np.linalg.norm(v-AA,axis=1)));allowed=normal[:,1]>=0;allowed[[start,end]]=True;use=allowed[edges].all(1);e=edges[use];length=np.linalg.norm(v[e[:,0]]-v[e[:,1]],axis=1);G=coo_matrix((np.r_[length,length],(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),shape=(len(v),len(v))).tocsr();distance,parent=dijkstra(G,indices=start,return_predecessors=True);assert np.isfinite(distance[end]),'posterior spine path unavailable';path=[end]
   while path[-1]!=start:path.append(int(parent[path[-1]]))
   curve=v[path]
  else:
   medial=np.linalg.norm(v-TS,axis=1)<=np.linalg.norm(v-AA,axis=1);superior=v[medial,2].max();curve=[]
   for z in np.arange(AI[2],superior+spacing,spacing):
    ids=np.flatnonzero(np.abs(v[:,2]-z)<=spacing)
    if len(ids):curve.append(v[ids[np.argmax(v[ids,0])]])
   curve=np.array(curve)
   if key.endswith('upper-medial'):curve=curve[curve[:,2]>=TS[2]]
   if key.endswith('lower-medial'):curve=curve[curve[:,2]<=TS[2]]
  distance=cKDTree(curve).query(v)[0];keep=(distance[base.faces]<=spacing).any(1);roi=trimesh.Trimesh(vertices=v.copy(),faces=base.faces[keep],process=False);roi.remove_unreferenced_vertices();assert len(roi.faces)>0;meshes[key]=roi
 if key not in meshes:
  file={'scapula':'FMA13395.stl','clavicle':'FMA13322.stl','humerus':'FMA23130.stl'}.get(key,byid.get(key));meshes[key]=trimesh.load_mesh(cache/file,process=True)
 return meshes[key]
originTargets={'v58RhomboidMajor':['T1','T2','T3','T4','T5'],'v58RhomboidMinor':['C7','T1'],'v59TrapMiddle':['C7','T1','T2','T3','T4'],'v60PecClav':['clavicle'],'v60DeltPost':['scapula']}
records=[];bindings=[]
for B in json.loads(src.read_text()):
 p=(np.array(B['rest']).reshape(-1,3)/scale-offset)@P;w=np.array(B['weights']);old=w.copy();f=np.array(B['faces']).reshape(-1,3);edges=np.unique(np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1),axis=0);spacing=float(np.median(np.linalg.norm(p[edges[:,0]]-p[edges[:,1]],axis=1)));ends=[]
 for end,target in enumerate([B['spec']['a'],B['spec']['b']]):
  roiTargets={'v58RhomboidMajor':['scapula:lower-medial'],'v58RhomboidMinor':['scapula:upper-medial'],'v58Levator':['scapula:upper-medial'],'v58Serratus':['scapula:medial'],'v59TrapMiddle':['scapula:spine'],'v59TrapLower':['scapula:spine']} if end==1 else{}
  if (w==end).any() and B['spec']['id']not in roiTargets:ends.append({'end':end,'status':'existing-region-retained','vertices':int((w==end).sum())});continue
  targets=originTargets.get(B['spec']['id'])if end==0 else roiTargets.get(B['spec']['id'],[target])
  if not targets:ends.append({'end':end,'status':'missing-anatomical-target-surface'});continue
  normalized=(old-old.min())/max(old.max()-old.min(),1e-20);candidate=np.flatnonzero(normalized<=.5 if end==0 else normalized>=.5);candidate=candidate[(old[candidate]!=1)if end==0 else(old[candidate]!=0)]
  if B['spec']['id']=='v58Serratus' and end==1:candidate=np.arange(len(w))
  else:candidate=candidate[w[candidate]!=(1-end)]
  distances=np.full(len(candidate),np.inf);nearest=np.zeros((len(candidate),3));targetIds=np.zeros(len(candidate),int);targetSpacing=[float(np.median(bone(key).edges_unique_length))for key in targets];resolution=spacing+float(np.median(targetSpacing));best=np.inf
  for j,key in enumerate(targets):
   m=bone(key);vertexDistance=cKDTree(m.vertices).query(p[candidate])[0];first=int(np.argmin(vertexDistance));_,firstDistance,_=trimesh.proximity.closest_point(m,p[candidate[first:first+1]]);best=min(best,float(firstDistance[0]));possible=np.flatnonzero(vertexDistance<=best+resolution+float(m.edges_unique_length.max()))
   # Vertex distance minus the longest source edge is a lower bound on
   # triangle distance. Prune only when a point cannot belong to this region.
   for start in range(0,len(possible),256):
    idx=possible[start:start+256];q,d,_=trimesh.proximity.closest_point(m,p[candidate[idx]]);take=d<distances[idx];selected=idx[take];distances[selected]=d[take];nearest[selected]=q[take];targetIds[selected]=j
  resolution=spacing+float(np.median(targetSpacing));minimum=float(distances.min());mask=distances<=minimum+resolution
  ids=candidate[mask]
  if len(ids)<3:ends.append({'end':end,'status':'insufficient-surface-region','vertices':len(ids)});continue
  overwritten=int((w[ids]==(1-end)).sum());w[ids]=end
  ends.append({'end':end,'status':'source-terminal-candidate','targetBones':targets,'vertices':len(ids),'overrodeOppositeLegacyPins':overwritten,'centroidAtlasMm':p[ids].mean(0).tolist(),'nearestBoneCentroidAtlasMm':nearest[mask].mean(0).tolist(),'minimumGapMm':minimum,'candidateResolutionMm':resolution,'maxGapMm':float(distances[mask].max()),'method':'bone-near terminal half-domain; tolerance from measured source mesh edge spacing','confidence':'estimated muscle/tendon terminal, not measured enthesis'})
 B['weights']=w.tolist();B['material']=w.tolist();B['terminalRegions']={'method':__doc__,'ends':ends};bindings.append(B);records.append({'id':B['spec']['id'],'ends':ends});print(json.dumps(records[-1]),flush=True)
(out/'bindings.json').write_text(json.dumps(bindings,separators=(',',':')));(out/'report.json').write_text(json.dumps({'method':__doc__,'muscles':records},indent=2)+'\n')
