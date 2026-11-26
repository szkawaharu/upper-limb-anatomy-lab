"""Ray-measured rest thickness on exact source surfaces; no invented shape.
Missing/ambiguous rays are explicitly marked and receive no thickness offset.
The ray half-chord is a geometric estimate, not histological belly thickness.
"""
import json,sys,hashlib
from pathlib import Path
import numpy as np,trimesh
src=Path(sys.argv[1]if len(sys.argv)>1 else'test-results/v076-source-bindings/bindings.json')
out=Path(sys.argv[2]if len(sys.argv)>2 else'test-results/v076-source-thickness');out.mkdir(parents=True,exist_ok=True)
records=[]
for B in json.loads(src.read_text()):
    p=np.array(B['rest']).reshape(-1,3);f=np.array(B['faces']).reshape(-1,3);scale=B['definition']['source']['worldUnitsPerMm'];p/=scale
    faceN=np.cross(p[f[:,1]]-p[f[:,0]],p[f[:,2]]-p[f[:,0]]);N=np.zeros(p.shape)
    for k in range(3):np.add.at(N,f[:,k],faceN)
    N/=np.maximum(np.linalg.norm(N,axis=1,keepdims=True),1e-20)
    mesh=trimesh.Trimesh(vertices=p,faces=f,process=False,vertex_normals=N)
    eps=np.linalg.norm(p.max(0)-p.min(0))*1e-6
    origins=p-N*eps
    hit,idx,_=mesh.ray.intersects_location(origins,-N,multiple_hits=False)
    distance=np.einsum('ij,ij->i',p[idx]-hit,N[idx]);valid=distance>2*eps
    h=np.zeros(len(p),dtype='<f4');h[idx[valid]]=distance[valid]*.5*scale
    raw=h.tobytes();(out/(B['spec']['id']+'.f32')).write_bytes(raw)
    r={'id':B['spec']['id'],'vertices':len(p),'sourceClosed':bool(mesh.is_watertight),'sourceOriented':bool(mesh.is_winding_consistent),'sourceVolumeMm3':float(mesh.volume),'missingOrAmbiguousRays':int((h==0).sum()),'halfChordMm':np.quantile(h[h>0]/scale,[0,.5,.95,1]).tolist(),'epsilonMm':float(eps),'sha256':hashlib.sha256(raw).hexdigest()};records.append(r);print(json.dumps(r),flush=True)
    (out/'report.json').write_text(json.dumps({'method':__doc__,'inputSha256':hashlib.sha256(src.read_bytes()).hexdigest(),'meshes':records},indent=2)+'\n')
