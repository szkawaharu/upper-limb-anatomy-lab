"""Read-only provenance for high-elevation local shape changes.
All element and surface ratios are relative to the exact source geometry.
No arbitrary physiological strain limit, no geometry edits or release claim."""
import json,os
from pathlib import Path
import numpy as np
from volumetric_reference import Solid
src=Path(os.environ['QA_INPUT']);out=Path(os.environ['QA_OUTPUT']);out.mkdir(parents=True,exist_ok=True)
files=sorted(src.glob('accepted-*-index-*.npz'));assert files
a=np.load(files[0]);rest=a['rest'];tets=a['tets'];faces=a['faces'];classes=a['classes'];center=a['center'];fixed=a['fixed']
solid=Solid(rest,tets,fixed)
edges=np.unique(np.sort(np.concatenate([faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]),axis=1),axis=0)
def area(p):
    v=p[faces];return np.linalg.norm(np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0]),axis=1)/2
def length(p):return np.linalg.norm(p[edges[:,1]]-p[edges[:,0]],axis=1)
area0=area(rest);edge0=length(rest);assert area0.min()>0 and edge0.min()>0
labels=['C1','C2','C3','C4','scapula']
# Every boundary face must be adjacent to exactly one source tetrahedron.
mapping={}
for k in range(4):
    for i,t in enumerate(np.sort(np.delete(tets,k,axis=1),axis=1)):
        key=tuple(t);mapping[key]=i if key not in mapping else -1
face_tet=np.array([mapping[tuple(sorted(f))] for f in faces]);assert np.all(face_tet>=0)
def distribution(v,w):
    return {'quantiles':np.quantile(v,[0,.001,.01,.5,.99,.999,1]).tolist(),'weightFractionOutsideHalfToDouble':float(w[(v<.5)|(v>2)].sum()/w.sum())}
rows=[]
for file in files:
    d=np.load(file)
    for k,v in [('rest',rest),('tets',tets),('faces',faces),('fixed',fixed),('classes',classes),('center',center)]:assert np.array_equal(d[k],v)
    z=d['nodes'];F=solid.edges(z)@solid.inverse;J=np.linalg.det(F);assert np.isfinite(J).all() and J.min()>1e-9
    stretches=np.linalg.svd(F,compute_uv=False);ar=area(z)/area0;er=length(z)/edge0
    def extrema(v,ids,reference,kind):
        result=[]
        for side,order in [('min',np.argsort(v)[:5]),('max',np.argsort(v)[-5:][::-1])]:
            for i in order:
                vertices=ids[i]
                row={'kind':kind,'extreme':side,'id':int(i),'ratio':float(v[i]),'vertexIds':vertices.tolist(),'fixedLabels':classes[vertices].tolist(),'sourceCentroidAtlasMm':(rest[vertices].mean(0)+center).tolist(),'referenceMeasure':float(reference[i]),'posedMeasure':float(reference[i]*v[i])}
                if kind=='face':
                    j=face_tet[i];row.update(tetId=int(j),tetJ=float(J[j]),tetPrincipalStretches=stretches[j].tolist(),tetReferenceVolumeMm3=float(solid.volumes[j]),tetReferenceCondition=float(np.linalg.cond(solid.Dm[j])))
                result.append(row)
        return result
    rows.append({'checkpoint':file.name,'index':int(d['index']),'area':distribution(ar,area0),'edge':distribution(er,edge0),'minimumPrincipalStretch':distribution(stretches[:,-1],solid.volumes),'maximumPrincipalStretch':distribution(stretches[:,0],solid.volumes),'maxConditionF':float(np.max(stretches[:,0]/stretches[:,-1])),'JRange':[float(J.min()),float(J.max())],'extrema':extrema(ar,faces,area0,'face')+extrema(er,edges,edge0,'edge')})
    row=rows[-1]
    print('LEVATOR_SHAPE_CASE '+json.dumps({'checkpoint':file.name,'index':row['index'],'area':row['area'],'edge':row['edge'],'maxConditionF':row['maxConditionF'],'extrema':[row['extrema'][0],row['extrema'][5],row['extrema'][10],row['extrema'][15]]}),flush=True)
report={'rows':rows,'referenceAreaMm2':float(area0.sum()),'referenceAreaQuantilesMm2':np.quantile(area0,[0,.01,.5,.99,1]).tolist(),'halfToDoubleIsDiagnosticOnly':True,'productAccepted':False}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print('LEVATOR_SHAPE_REPORT '+json.dumps({'poses':len(rows),'maxConditionF':max(r['maxConditionF'] for r in rows),'maxAreaFractionOutsideHalfToDouble':max(r['area']['weightFractionOutsideHalfToDouble'] for r in rows)}),flush=True)
