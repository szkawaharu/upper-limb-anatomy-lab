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
