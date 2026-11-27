"""Experimental local log-volume constraints on the audited intersection cage.

Not a production engine or anatomical acceptance test. No cell is discarded
by size. Background integration and original fixed nodes are unchanged.
A P1 displacement / independent P0 cell constraint can be overdetermined;
pressure continuation alone is not a remedy for volumetric locking.
"""
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.csgraph import structural_rank
from scipy.sparse.linalg import lsqr


def log_volume_jacobian(solid, x, cells):
    """Exact log(det F) and derivative in full node-major coordinates."""
    cells = np.asarray(cells, dtype=int)
    F = solid.edges(x)[cells] @ solid.inverse[cells]
    J = np.linalg.det(F)
    if not np.isfinite(J).all() or np.any(J <= 1e-9):
        raise ValueError('Non-positive local volume constraint input')
    inverse = solid.inverse[cells]
    gradients = np.concatenate([-inverse.sum(1)[:, None, :], inverse], axis=1)
    values = np.einsum('tij,taj->tai', np.linalg.inv(F).transpose(0, 2, 1), gradients)
    columns = (solid.tets[cells, :, None]*3 + np.arange(3)).reshape(-1, 12)
    rows = np.broadcast_to(np.arange(len(cells))[:, None], columns.shape)
    A = coo_matrix((values.ravel(), (rows.ravel(), columns.ravel())),
                   shape=(len(cells), x.size)).tocsr()
    A.eliminate_zeros()
    return np.log(J), A


def physical_cells(solid, source_volumes, background_weights):
    source = np.asarray(source_volumes, float)
    background = np.asarray(background_weights, float)
    if (source.shape != solid.weights.shape or background.shape != source.shape
            or not np.isfinite(source).all() or not np.isfinite(background).all()
            or np.any(source < 0) or np.any(background < 0)
            or source.sum() <= 0 or background.sum() >= 1):
        raise ValueError('Invalid audited material/background integration')
    w = source/source.sum()*(1-background.sum())
    if not np.allclose(w+background, solid.weights, rtol=1e-10, atol=1e-14):
        raise ValueError('Material/background integration does not reconstruct weights')
    cells = np.flatnonzero(source > 0)
    return cells, w[cells]


def audit(solid, x, cells, max_iterations=2000):
    """Unweighted linearized feasibility, not a nonlinear impossibility proof.
    Structural rank is only an upper bound on numerical rank. LSQR residual is
    from a trial least-squares direction; iteration-limit is explicitly stored.
    """
    c, A = log_volume_jacobian(solid, x, cells)
    free = (solid.free[:, None]*3+np.arange(3)).ravel()
    B = A[:, free].tocsr()
    B.eliminate_zeros()
    # A zero differential at one configuration is not a proof that a
    # nonlinear constraint is immovable. Only cells with every node fixed
    # are strictly immutable under the stated Dirichlet data.
    zero_rows = np.diff(B.indptr) == 0
    immovable = np.all(~np.isin(solid.tets[np.asarray(cells)], solid.free), axis=1)
    if len(free):
        # Column equilibration changes the search coordinates, not the
        # unweighted residual or the linearized feasible subspace.
        norms = np.sqrt(np.asarray(B.power(2).sum(axis=0)).ravel())
        scale = np.ones_like(norms)
        np.divide(1., norms, out=scale, where=norms>0)
        fit = lsqr(B@diags(scale), -c, atol=1e-11, btol=1e-11, iter_lim=max_iterations)
        residual = B@(scale*fit[0])+c
        stop, count = int(fit[1]), int(fit[2])
        stationarity = float(np.linalg.norm(B.T@residual))
    else:
        residual = c
        stop, count, stationarity = 0, 0, 0.
    rank_bound = int(structural_rank(B)) if len(free) else 0
    return {'physicalCells':len(c), 'freeDisplacementDofs':len(free),
            'structuralRankUpperBound':rank_bound,
            'constraintRankDeficiencyLowerBound':len(c)-rank_bound,
            'immovableCellIds':np.asarray(cells)[immovable].tolist(),
            'zeroLinearizedRowCount':int(zero_rows.sum()),
            'zeroLinearizedButNotFixedCount':int((zero_rows & ~immovable).sum()),
            'lsqrColumnEquilibrated':True,
            'lsqrIterationLimit':max_iterations,
            'immovableMaxAbsLogJ':float(np.max(np.abs(c[immovable]), initial=0)),
            'maxAbsLogJ':float(np.max(np.abs(c), initial=0)),
            'linearizedResidualL2':float(np.linalg.norm(residual)),
            'linearizedResidualMax':float(np.max(np.abs(residual), initial=0)),
            'linearizedNormalResidualL2':stationarity,
            'lsqrStop':stop, 'lsqrIterations':count,
            'nonlinearFeasibilityEstablished':False}


class LocalPressure:
    """Augmented Lagrangian term on physical log-volume only.
    Exact scalar objective and gradient. PSD Gauss-Newton search metric omits
    log-volume curvature; the unchanged positive-J Armijo guard decides steps.
    """
    def __init__(self, solid, cells, weights, pressure, rho):
        self.solid, self.cells = solid, np.asarray(cells)
        self.weights, self.pressure = np.asarray(weights), np.asarray(pressure).copy()
        self.rho = float(rho)
        if (self.pressure.shape != self.weights.shape or self.cells.shape != self.weights.shape
                or not np.isfinite(self.pressure).all() or not np.isfinite(self.rho)
                or self.rho <= 0 or not np.all(self.weights > 0)):
            raise ValueError('Invalid local pressure input')

    def evaluate(self, x, gradient=True):
        c, A = log_volume_jacobian(self.solid, x, self.cells)
        energy = np.dot(self.weights, self.pressure*c + .5*self.rho*c*c)
        g = (A.T@(self.weights*(self.pressure+self.rho*c))).reshape(x.shape) if gradient else None
        return float(energy), g

    def hessian(self, x):
        _, A = log_volume_jacobian(self.solid, x, self.cells)
        return A.T@diags(self.rho*self.weights)@A


def equilibrate(solid, initial, contact, target, source_volumes, background_weights,
                *, max_log_error, rho=30., bulk=30., iterations=150, outer_steps=20):
    """Experimental fixed-pose solve; no relaxed geometry or fallback acceptance.
    Caller supplies an explicit local tolerance for the experiment, not a new
    product gate. A failed outer solve returns the original accepted input.
    No pressure or partial solution leaks into the next joint pose.
    """
    from actual_bone_contact import ActualDistancePotential
    from contact_potential import AugmentedSolid
    from cage_solver import minimize
    if not np.isfinite(max_log_error) or max_log_error <= 0 or outer_steps < 1:
        raise ValueError('An explicit positive local experimental tolerance is required')
    cells, weights = physical_cells(solid, source_volumes, background_weights)
    x = initial.copy()
    pins = x[solid.fixed].copy()
    pressure = np.zeros(len(cells))
    penalty = contact.initial_penalty
    rows = []
    before = audit(solid, x, cells)
    def failure(reason):
        return initial.copy(), {'converged':False, 'accepted':False, 'reason':reason,
                                'initialAudit':before, 'history':rows}
    if before['immovableMaxAbsLogJ'] > max_log_error:
        return failure('fixed-local-volume-conflict')
    for outer in range(outer_steps):
        local = LocalPressure(solid, cells, weights, pressure, rho)
        bone = ActualDistancePotential(contact, len(x), target, penalty)
        combined = AugmentedSolid(solid, [local, bone])
        trial, report = minimize(combined, x, bulk=bulk, iterations=iterations)
        if not np.array_equal(trial[solid.fixed], pins):
            raise RuntimeError('Local pressure solve changed fixed positions')
        _, gradient, J = combined.evaluate(trial, bulk)
        residual = float(np.linalg.norm(solid.factor.solve(gradient[solid.free]), axis=1).max(initial=0))
        c = np.log(J[cells])
        maximum = float(np.max(np.abs(c)))
        depth = max(r['maxDepthIncreaseMm'] for r in bone.state(trial)[1])
        rows.append({'outer':outer, 'maxAbsLogJ':maximum,
                     'physicalJRange':[float(J[cells].min()),float(J[cells].max())],
                     'maxDepthIncreaseMm':depth, 'residualMm':residual,
                     'contactPenalty':penalty, 'pressureMaxAbs':float(np.max(np.abs(pressure))),
                     'innerConverged':report['converged']})
        if not report['converged'] or residual >= .001 or not np.isfinite(J).all() or J.min() <= 1e-9:
            return failure('local-pressure-inner-failed')
        if maximum <= max_log_error and depth <= .1:
            return trial, {'converged':True, 'accepted':False,
                           'reason':'experimental-local-volume-and-contact',
                           'initialAudit':before, 'finalAudit':audit(solid, trial, cells),
                           'history':rows, 'experimentalMaxAbsLogJ':max_log_error}
        x = trial
        pressure += rho*c
        if depth > .1:
            penalty *= 10
    return failure('local-pressure-outer-limit')


if __name__ == '__main__':
    import argparse, hashlib, json
    from pathlib import Path
    from volumetric_reference import Solid
    parser = argparse.ArgumentParser(description='Read-only local constraint audit of an integrated cage and saved pose')
    parser.add_argument('cage', type=Path)
    parser.add_argument('snapshot', type=Path)
    args = parser.parse_args()
    with np.load(args.cage) as a, np.load(args.snapshot) as s:
        if not (np.array_equal(a['nodes'], s['referenceNodes']) and np.array_equal(a['tets'], s['cageTets'])):
            raise ValueError('Snapshot and cage geometry differ')
        solid = Solid(a['nodes'], a['tets'], np.flatnonzero(a['classes'] >= 0), a['weights'])
        cells, _ = physical_cells(solid, a['sourceMaterialVolumes'], a['backgroundWeights'])
        result = audit(solid, s['nodes'], cells)
    result['inputsSha256'] = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (args.cage,args.snapshot)}
    result.update(experimental=True, accepted=False)
    print(json.dumps(result, indent=2))
