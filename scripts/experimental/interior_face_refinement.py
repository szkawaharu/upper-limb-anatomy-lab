"""Conforming interior-face enrichment for a full physical P1 solid.

A tetrahedron-centre bubble alone cannot change its integrated signed volume
if its boundary is fixed. A shared interior-face centre provides a normal
flux between its two incident parents. This is a targeted h-refinement
experiment, not an inf-sup-stable mixed FEM or a production model.

Only identity/full-source material integration is supported. A cut-cell cage
would require new audited intersections rather than inherited densities.
"""
from itertools import combinations
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from volumetric_reference import Solid


def face_owners(tets):
    owners = {}
    for i, tet in enumerate(tets):
        for face in combinations(map(int,tet),3):
            owners.setdefault(tuple(sorted(face)),[]).append(i)
    if any(len(v)>2 for v in owners.values()):
        raise ValueError('Nonmanifold tetrahedral face')
    return owners


def select_faces(solid, x, budget):
    """Rank eligible interior faces by maximum incident parent abs(log J).
    No surface/fixed face is split; no fully fixed parent is refined.
    The budget is a declared experiment size, not an exclusion from J audit.
    """
    if budget<1: raise ValueError('Positive face budget required')
    J=np.linalg.det(solid.edges(x)@solid.inverse)
    if np.any(J<=0) or not np.isfinite(J).all(): raise ValueError('Invalid pose')
    score=np.abs(np.log(J))
    fixed=set(map(int,solid.fixed))
    full=np.isin(solid.tets,solid.fixed).all(axis=1)
    owners=face_owners(solid.tets)
    eligible=[(f,ids) for f,ids in owners.items() if len(ids)==2 and
              not any(full[ids]) and not all(i in fixed for i in f)]
    eligible.sort(key=lambda item:(-max(score[item[1]]),item[0]))
    return [f for f,_ in eligible[:budget]]


def refine(solid, x, faces):
    """Add shared face centres and affected parent centres, with exact
    piecewise-affine transfer. External triangles and pins stay identical.
    All newly created nodes are internal to the original solid.
    """
    if not np.allclose(solid.weights,solid.volumes/solid.volumes.sum(),rtol=1e-11,atol=1e-15):
        raise ValueError('Cut-cell/nonuniform material requires independent integration')
    owners=face_owners(solid.tets)
    faces=sorted(set(tuple(sorted(map(int,f))) for f in faces))
    fixed=set(map(int,solid.fixed))
    full=np.isin(solid.tets,solid.fixed).all(axis=1)
    if not faces: raise ValueError('No interior faces selected')
    if any(len(owners.get(f,[]))!=2 or all(i in fixed for i in f)
           or any(full[owners[f]]) for f in faces):
        raise ValueError('Exterior or fully fixed face/parent cannot be enriched')
    rest=list(solid.rest.copy()); pose=list(np.asarray(x).copy())
    def add(ids):
        idx=len(rest);rest.append(solid.rest[list(ids)].mean(0));pose.append(x[list(ids)].mean(0));return idx
    fnode={f:add(f) for f in faces}
    affected=sorted({i for f in faces for i in owners[f]})
    cnode={i:add(solid.tets[i]) for i in affected}
    children=[];parent=[]
    for i,tet in enumerate(solid.tets):
        if i not in cnode:
            children.append(tuple(tet));parent.append(i);continue
        for f in combinations(map(int,tet),3):
            f=tuple(sorted(f))
            triangles=[f] if f not in fnode else [(f[j],f[(j+1)%3],fnode[f]) for j in range(3)]
            for tri in triangles:
                children.append((cnode[i],*tri));parent.append(i)
    rest=np.asarray(rest);pose=np.asarray(pose);tets=np.asarray(children,dtype=np.int32);parent=np.asarray(parent)
    p=rest[tets];det=np.linalg.det(np.stack([p[:,j]-p[:,0] for j in (1,2,3)],axis=2))
    flip=det<0;tets[flip,2],tets[flip,3]=tets[flip,3].copy(),tets[flip,2].copy()
    refined=Solid(rest,tets,solid.fixed)
    covered=np.bincount(parent,weights=refined.volumes,minlength=len(solid.tets))
    if np.max(np.abs(covered/solid.volumes-1))>1e-8:
        raise ValueError('Parent volume coverage mismatch')
    external={f for f,v in owners.items() if len(v)==1}
    after={f for f,v in face_owners(tets).items() if len(v)==1}
    if external!=after: raise ValueError('Refinement changed exterior triangulation')
    assert np.array_equal(refined.rest[:len(solid.rest)],solid.rest)
    assert np.array_equal(refined.fixed,solid.fixed)
    assert np.array_equal(pose[:len(x)],x)
    # Signed parent-volume flux via normal motion of the shared face centre.
    # Independent incidence columns sum to zero; a tree gives n-1 directions
    # per component. This does not certify numerical rank of the combined
    # original-plus-enriched Jacobian.
    rows=[];cols=[];values=[]
    for j,f in enumerate(faces):
        a,b=owners[f];rows.extend([a,b]);cols.extend([j,j]);values.extend([1,-1])
    incidence=coo_matrix((values,(rows,cols)),shape=(len(solid.tets),len(faces))).tocsr()
    sub=incidence[affected];components,_=connected_components(sub@sub.T,directed=False)
    info={'selectedInteriorFaces':len(faces),'affectedParents':len(affected),
          'affectedParentIds':affected,'addedNodes':len(rest)-len(solid.rest),
          'newTetrahedra':len(tets),'maxRelativeParentCoverageError':float(np.max(np.abs(covered/solid.volumes-1))),
          'exteriorFacesUnchanged':True,'fixedNodesUnchanged':True,
          'incidenceFluxRank':len(affected)-components,
          'incidenceRankIsNotCombinedJacobianRank':True}
    return refined,pose,parent,info


def metrics(original, refined, x, parent):
    J=np.linalg.det(refined.edges(x)@refined.inverse)
    ratio=np.bincount(parent,weights=refined.volumes*J,minlength=len(original.tets))/original.volumes
    if not np.isfinite(J).all() or np.any(J<=1e-9):
        raise ValueError('Nonpositive refined Jacobian')
    return {'parentJRange':[float(ratio.min()),float(ratio.max())],
            'allChildJRange':[float(J.min()),float(J.max())],
            'maxParentAbsLogJ':float(np.max(np.abs(np.log(ratio)))),
            'maxChildAbsLogJ':float(np.max(np.abs(np.log(J)))),
            'globalVolumeRatio':float(np.dot(refined.weights,J)),
            'parentJ':ratio}
