"""Source-derived continuous Dirichlet basis for a reduced solid model.
Control patches are numerical DOFs, not anatomical attachments. Harmonic
functions vanish exactly on the original fixed regions; affine enrichment
allows patches to rotate/shear, rather than making them rigid translations.
"""
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.csgraph import connected_components,dijkstra
from scipy.sparse.linalg import splu


def build_basis(nodes,tets,source_weights,controls=24):
    n=len(nodes);source_count=len(source_weights);active=np.unique(tets)
    fixed=np.flatnonzero((source_weights==0)|(source_weights==1))
    free=np.setdiff1d(active,fixed)
    edges=np.unique(np.sort(np.concatenate([tets[:,[i,j]]for i in range(4)for j in range(i+1,4)]),axis=1),axis=0)
    length=np.linalg.norm(nodes[edges[:,0]]-nodes[edges[:,1]],axis=1)
    if np.any(length<=0):raise ValueError('Zero source tetrahedral edge')
    a=np.r_[edges[:,0],edges[:,1]];b=np.r_[edges[:,1],edges[:,0]];v=np.r_[length,length]
    distances=coo_matrix((v,(a,b)),shape=(n,n)).tocsr()
    adjacency=coo_matrix((1/v,(a,b)),shape=(n,n)).tocsr();L=diags(np.asarray(adjacency.sum(1)).ravel())-adjacency
    _,labels=connected_components(adjacency)
    for k in np.unique(labels[free]):
        if not np.any(labels[fixed]==k):raise ValueError('Unanchored source component; attachment data missing')
    count=min(controls,len(free));chosen=[];distance_fields=[];nearest=np.full(n,np.inf)
    current=int(free[np.argmin(np.linalg.norm(nodes[free]-nodes[free].mean(0),axis=1))])
    for _ in range(count):
        chosen.append(current);field=dijkstra(distances,indices=current);distance_fields.append(field);nearest=np.minimum(nearest,field);nearest[chosen]=-1
        current=int(free[np.argmax(nearest[free])])
    chosen=np.array(chosen,dtype=int);control_distance=np.array(distance_fields)
    pin_distance=dijkstra(distances,indices=fixed,min_only=True)
    pair=control_distance[:,chosen].copy();np.fill_diagonal(pair,np.inf)
    radius=.5*np.minimum(pair.min(1),pin_distance[chosen])
    if not np.all(np.isfinite(radius)&(radius>0)):raise ValueError('Undefined control patch extent')
    normalized=control_distance/radius[:,None];owner=np.argmin(normalized,axis=0);inside=normalized[owner,np.arange(n)]<=1+1e-12
    inside[fixed]=False;inside[np.setdiff1d(np.arange(n),active)]=False
    patch=np.flatnonzero(inside);unknown=np.setdiff1d(active,np.r_[fixed,patch])
    phi=np.zeros((n,count));phi[patch,owner[patch]]=1
    if len(unknown):
        factor=splu(L[unknown][:,unknown].tocsc());phi[unknown]=factor.solve(-(L[unknown][:,patch]@phi[patch]))
    if phi.min()<-1e-9 or phi.max()>1+1e-9 or phi.sum(1).max()>1+1e-9:raise ValueError('Harmonic maximum principle violated')
    # A surface-only harmonic field can disagree across the two faces of a
    # thin solid. Preserve the actual fixed regions, and solve the NUMERICAL
    # lifting coordinate on the full volume graph. This is not a fascicle arc
    # coordinate or an active contraction value. Attachment masks stay separate.
    weights=np.zeros(n);weights[fixed]=source_weights[fixed]
    if len(free):weights[free]=splu(L[free][:,free].tocsc()).solve(-(L[free][:,fixed]@weights[fixed]))
    if weights[active].min()<-1e-9 or weights[active].max()>1+1e-9:raise ValueError('Volumetric lifting coordinate outside its Dirichlet range')
    raw=np.zeros((n,count*4))
    for i,node in enumerate(chosen):
        raw[:,i*4]=phi[:,i]
        raw[:,i*4+1:i*4+4]=phi[:,i,None]*(nodes-nodes[node])/radius[i]
    p=nodes[tets];Dm=np.stack([p[:,j]-p[:,0]for j in[1,2,3]],axis=2);volume=np.linalg.det(Dm)/6
    if np.any(volume<=0):raise ValueError('Non-positive reference volume')
    inverse=np.linalg.inv(Dm);grad=np.concatenate([-inverse.sum(1)[:,None,:],inverse],axis=1)
    G=np.einsum('tij,tik->tjk',grad,raw[tets]);weighted=(G*np.sqrt(volume/volume.sum())[:,None,None]).reshape(-1,raw.shape[1]);S=weighted.T@weighted
    eigen,U=np.linalg.eigh(S);keep=eigen>eigen.max()*1e-12
    transform=U[:,keep]/np.sqrt(eigen[keep]);basis=raw@transform
    if np.any(basis[fixed]!=0):raise ValueError('Attachment boundary changed')
    # A scalar skin weight times x is not a harmonic vector extension on a
    # sliver tetrahedron. Extend the FOUR affine boundary channels together in
    # the actual FEM metric. Dynamic lifting is A(x)+H_B[(B-A)(x)], so common
    # rigid motion is exact and both wide fixed regions remain exact.
    fem_local=np.einsum('tij,tkj,t->tik',grad,grad,volume/volume.sum())
    rows=np.broadcast_to(tets[:,:,None],fem_local.shape).ravel();cols=np.broadcast_to(tets[:,None,:],fem_local.shape).ravel()
    fem=coo_matrix((fem_local.ravel(),(rows,cols)),shape=(n,n)).tocsr()
    lifting=np.zeros((n,4));insertion=fixed[source_weights[fixed]==1];lifting[insertion,:3]=nodes[insertion];lifting[insertion,3]=1
    if len(free):lifting[free]=splu(fem[free][:,free].tocsc()).solve(-(fem[free][:,fixed]@lifting[fixed]))
    return {'basis':basis,'weights':weights,'lifting':lifting,'liftingMethod':np.array('fem-dirichlet-affine-boundary'),'originalSurfaceWeights':source_weights,'fixed':fixed,'controls':chosen,'radiusMm':radius,'controlPatchVertices':np.bincount(owner[patch],minlength=count),'eigenvalues':eigen,'discardedNumericalModes':int((~keep).sum()),'referenceInverse':inverse,'referenceVolume':volume,'basisTransform':transform}
