"""Real-source, fixed-bone-boundary local-volume experiment helpers.
No mesh contact oracle is supplied in this fixture. These experiments are
explicitly material-only, not acceptance and not a deltoid intersection run.
"""
import numpy as np
from local_volume_constraints import LocalPressure, audit, physical_cells
from contact_potential import AugmentedSolid
from cage_solver import minimize


def inspect_pose(solid, x, iterations=128):
    cells, weights = physical_cells(solid, solid.volumes, np.zeros_like(solid.weights))
    diagnostic = audit(solid, x, cells, max_iterations=iterations)
    # Store full fixed IDs in the artifact only via a compact count and
    # representative worst cell; the counts cover every original tetrahedron.
    ids = np.asarray(diagnostic.pop('immovableCellIds'), dtype=int)
    J = np.linalg.det(solid.edges(x) @ solid.inverse)
    diagnostic.update(fullyFixedCells=int(len(ids)),
                      fixedWorstCellId=int(ids[np.argmax(np.abs(np.log(J[ids])))]) if len(ids) else None,
                      JRange=[float(J.min()),float(J.max())],
                      volumeWeightedJ=float(np.dot(solid.weights,J)),
                      linearizedOnly=True, contactEvaluated=False)
    return diagnostic


def compare_pressure(solid, initial, *, rho=30., bulk=30.):
    """Equal first energy, then compare pressure update vs a larger bulk.
    Both branches start at the exact same converged first-stage geometry.
    Failed trials are recorded, never used as the next initial condition.
    """
    cells, weights = physical_cells(solid, solid.volumes, np.zeros_like(solid.weights))
    pins = initial[solid.fixed].copy()
    zero = np.zeros(len(cells))
    result = {'experimental':True,'accepted':False,'boneContactEvaluated':False,
              'selfCollisionEvaluated':False,'baseBulk':bulk,'rho':rho,
              'initial':inspect_pose(solid, initial),'stages':[]}
    baseline = AugmentedSolid(solid,[LocalPressure(solid,cells,weights,zero,rho)])
    eb, gb, _ = baseline.evaluate(initial,bulk)
    ec, gc, _ = solid.evaluate(initial,bulk+rho)
    assert np.isclose(eb,ec,rtol=1e-11,atol=1e-13)
    assert np.allclose(gb,gc,rtol=1e-10,atol=1e-12)
    result['zeroPressureMatchesBulkEnergyGradient']=True

    def run(label, model, start, coefficient):
        x, report = minimize(model, start, bulk=coefficient, iterations=80)
        assert np.array_equal(x[solid.fixed],pins)
        J = np.linalg.det(solid.edges(x) @ solid.inverse)
        assert np.isfinite(J).all() and J.min()>1e-9
        row = {'label':label,'converged':report['converged'],'reason':report['reason'],
               'iterations':report['iterations'],'residualMm':report['residualMaxMm'],
               'JRange':[float(J.min()),float(J.max())],
               'maxAbsLogJ':float(np.max(np.abs(np.log(J)))),
               'volumeWeightedJ':float(np.dot(solid.weights,J)),
               'maxDisplacementFromInitialMm':float(np.linalg.norm(x-initial,axis=1).max()),
               'attachmentMaxMm':float(np.linalg.norm(x[solid.fixed]-pins,axis=1).max(initial=0)),
               'usedAsNextInitial':False}
        result['stages'].append(row)
        return x,row

    first,row=run('equal-first-bulk',solid,initial,bulk+rho)
    if not row['converged'] or row['residualMm']>=.001:
        result['stopReason']='first-control-not-converged'
        return result
    row['usedAsNextInitial']=True
    pressure = rho*np.log(np.linalg.det(solid.edges(first)@solid.inverse))
    augmented=AugmentedSolid(solid,[LocalPressure(solid,cells,weights,pressure,rho)])
    _,p=run('local-pressure-update',augmented,first,bulk)
    _,c=run('common-bulk-control',solid,first,bulk+2*rho)
    result['comparisonConverged']=bool(p['converged'] and c['converged']
                                      and p['residualMm']<.001 and c['residualMm']<.001)
    if result['comparisonConverged']:
        result['pressureMinusControlMaxAbsLogJ']=p['maxAbsLogJ']-c['maxAbsLogJ']
    return result


def compare_refinement(solid, initial, surface_faces, budget=128, bulk=30.):
    from interior_face_refinement import refine, select_faces, metrics
    faces=select_faces(solid,initial,budget)
    refined,y,parent,info=refine(solid,initial,faces)
    initial_metrics=metrics(solid,refined,y,parent)
    initial_parent=initial_metrics.pop('parentJ')
    original_J=np.linalg.det(solid.edges(initial)@solid.inverse)
    assert np.allclose(initial_parent,original_J,rtol=1e-8,atol=1e-10)
    assert np.isclose(refined.evaluate(y,bulk,False)[0],solid.evaluate(initial,bulk,False)[0],rtol=1e-8,atol=1e-12)
    trial,report=minimize(refined,y,bulk=bulk,iterations=80)
    assert np.array_equal(trial[solid.fixed],initial[solid.fixed])
    final_metrics=metrics(solid,refined,trial,parent)
    final_parent=final_metrics.pop('parentJ')
    ids=np.asarray(info['affectedParentIds'])
    mask=np.ones(len(original_J),dtype=bool);mask[ids]=False
    volume=lambda v:float(np.einsum('ti,ti->t',v[surface_faces[:,0]],np.cross(v[surface_faces[:,1]],v[surface_faces[:,2]])).sum()/6)
    surface_ratio=volume(trial)/volume(solid.rest)
    assert abs(surface_ratio-final_metrics['globalVolumeRatio'])<1e-8
    return {'experimental':True,'accepted':False,'boneContactEvaluated':False,'selfCollisionEvaluated':False,
            'bulk':bulk,'faceBudget':budget,'refinement':info,'before':initial_metrics,'after':final_metrics,
            'converged':report['converged'],'residualMm':report['residualMaxMm'],'reason':report['reason'],
            'initialEnergy':solid.evaluate(initial,bulk,False)[0],'finalEnergy':report['energy'],
            'attachmentMaxMm':0.,'surfaceVolumeRatio':surface_ratio,
            'affectedParentBeforeJRange':[float(original_J[ids].min()),float(original_J[ids].max())],
            'affectedParentAfterJRange':[float(final_parent[ids].min()),float(final_parent[ids].max())],
            'unaffectedParentAfterJRange':[float(final_parent[mask].min()),float(final_parent[mask].max())] if mask.any() else None,
            'maxOriginalNodeDisplacementMm':float(np.linalg.norm(trial[:len(initial)]-initial,axis=1).max()),
            'comparisonOnly':True,'noSnapshotPromoted':True}


def compare_mixed_pressure(solid, initial):
    from mixed_pressure_reference import MixedPressure
    # Same original source mesh, fixed positions and local shear model.
    # Zero-jump pressure is identical to bulk90; compare the equilibria
    # with a single declared dimensionless stabilization, not a fitted muscle parameter.
    control=MixedPressure(solid,bulk_pressure=60.,jump=0.)
    mixed=MixedPressure(solid,bulk_pressure=60.,jump=1.)
    x,base=minimize(solid,initial,bulk=90.,iterations=80)
    assert np.array_equal(x[solid.fixed],initial[solid.fixed])
    result={'experimental':True,'accepted':False,'boneContactEvaluated':False,
            'selfCollisionEvaluated':False,'initial':control.audit(initial),
            'control':{'converged':base['converged'],'residualMm':base['residualMaxMm'],
                       'audit':control.audit(x)},'comparisonConverged':False}
    if not base['converged'] or base['residualMaxMm']>=.001:
        result['reason']='control-not-converged'
        return result
    y,report=mixed.minimize(x,bulk=30.,iterations=40)
    assert np.array_equal(y[solid.fixed],initial[solid.fixed])
    result['mixed']={**report,'audit':mixed.audit(y),
                     'attachmentMaxMm':float(np.linalg.norm(y[solid.fixed]-initial[solid.fixed],axis=1).max(initial=0))}
    result['comparisonConverged']=bool(report['converged'] and report['residualMm']<.001)
    if result['comparisonConverged']:
        a=result['control']['audit'];b=result['mixed']['audit']
        result['mixedMinusControlMaxAbsLogJ']=b['maxAbsLogJ']-a['maxAbsLogJ']
        result['mixedMinusControlPressureJumpRms']=b['pressureJumpRms']-a['pressureJumpRms']
    return result
