"""Transfer a source-conforming solution between internal tetrahedralizations.
Anatomical boundary vertices are copied exactly. Interior nodes use actual
containing old elements, never nearest-surface snapping or invented offsets.
"""
import sys,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

def transfer(reference,cells,current,new_reference,source_count,new_cells):
    if not np.array_equal(new_reference[:source_count],reference[:source_count]):
        raise ValueError('Source boundary coordinates/order changed')
    old=reference[cells];inverse=np.linalg.inv(np.stack([old[:,j]-old[:,0]for j in[1,2,3]],axis=2));tree=cKDTree(old.mean(1));out=new_reference.copy();out[:source_count]=current[:source_count]
    remaining=np.setdiff1d(np.unique(new_cells),np.arange(source_count));count=len(remaining);max_error=0.;full_search=0
    for k in [16,64,256,1024,len(cells)]:
        if not len(remaining):break
        k=min(k,len(cells));next_ids=[];batch=1 if k==len(cells)else 128
        for begin in range(0,len(remaining),batch):
            ids=remaining[begin:begin+batch];points=new_reference[ids];_,nearest=tree.query(points,k=k);nearest=np.asarray(nearest).reshape(len(ids),k);d=points[:,None,:]-old[nearest,0];v=np.einsum('nkij,nkj->nki',inverse[nearest],d);weights=np.concatenate([1-v.sum(2,keepdims=True),v],axis=2);inside=(weights.min(2)>=-1e-8)&(weights.max(2)<=1+1e-8)
            for j,id in enumerate(ids):
                valid=np.flatnonzero(inside[j])
                if not len(valid):next_ids.append(id);continue
                index=valid[0];tet=nearest[j,index];w=weights[j,index];out[id]=w@current[cells[tet]];error=np.linalg.norm(w@reference[cells[tet]]-new_reference[id]);max_error=max(max_error,float(error));full_search+=int(k==len(cells))
        remaining=np.array(next_ids,dtype=int)
    if len(remaining):raise ValueError('Interior node not in original source solid: '+repr(remaining.tolist()))
    return out,{'transferredInteriorNodes':count,'fullSearches':full_search,'maxReferenceReconstructionWorld':max_error,'sourceBoundaryDeltaWorld':float(np.linalg.norm(out[:source_count]-current[:source_count],axis=1).max())}

if __name__=='__main__':
    oldfile=Path(sys.argv[1]);newfile=Path(sys.argv[2]);sourcefile=Path(sys.argv[3]);output=Path(sys.argv[4]);old=np.load(oldfile);new=np.load(newfile);source=np.load(sourcefile);rest=new['nodes'];cells=new['tets'];current,report=transfer(old['referenceWorld'],old['tets'],old['world'],rest,len(source['weights']),cells)
    edges=lambda x:np.stack([x[cells[:,j]]-x[cells[:,0]]for j in[1,2,3]],axis=2)
    J=np.linalg.det(edges(current))/np.linalg.det(edges(rest));report.update({'minJ':float(J.min()),'maxJ':float(J.max()),'nonPositiveElements':int((J<=0).sum()),'old':str(oldfile),'new':str(newfile),'accepted':False})
    output.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    if np.any(J<=0):raise ValueError('Changed internal discretization produced non-positive initial elements')
    np.savez(output.with_suffix('.npz'),world=current,tets=cells,referenceWorld=rest,sourceFaces=source['faces'],fixed=old['fixed']);print(json.dumps(report),flush=True)
