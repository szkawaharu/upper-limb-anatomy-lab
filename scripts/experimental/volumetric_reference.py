"""Offline geometric continuum reference, not a real-time/clinical model.

Original surface nodes are retained. The energy is normalized isochoric
distortion plus a log-volume penalty; coefficients have no force units and
are not measured muscle material constants. Every line-search trial rejects
non-positive element Jacobians. Fixed source regions are exact Dirichlet data.
No TetGen library is linked here; the input is an audited generated NPZ mesh.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu


class Solid:
    def __init__(self, nodes, tets, fixed):
        self.rest = np.array(nodes, float)
        self.tets = np.array(tets, np.int32)
        self.fixed = np.array(fixed, np.int32)
        active = np.unique(self.tets)
        self.inactive = np.setdiff1d(np.arange(len(nodes)), active)
        self.free = np.setdiff1d(active, self.fixed)
        self.Dm = self.edges(self.rest)
        determinant = np.linalg.det(self.Dm)
        if not np.all(determinant > 0):
            raise ValueError('Reference tetrahedra must have positive signed volume')
        self.inverse = np.linalg.inv(self.Dm)
        self.volumes = determinant / 6
        self.weights = self.volumes / self.volumes.sum()
        # Scalar FEM stiffness is a source-derived, fixed SPD preconditioner.
        gradients = np.concatenate([-self.inverse.sum(1)[:, None, :], self.inverse], axis=1)
        local = np.einsum('tij,tkj,t->tik', gradients, gradients, self.weights)
        rows = np.broadcast_to(self.tets[:, :, None], local.shape).ravel()
        cols = np.broadcast_to(self.tets[:, None, :], local.shape).ravel()
        L = coo_matrix((local.ravel(), (rows, cols)), shape=(len(nodes), len(nodes))).tocsr()
        count, labels = connected_components(L)
        missing = [k for k in np.unique(labels[active]) if not np.any(labels[self.fixed] == k)]
        if missing:
            raise ValueError('Components without attachment data: ' + repr(missing))
        self.factor = splu(L[self.free][:, self.free].tocsc())
        self.coupling = L[self.free][:, self.fixed]

    def edges(self, x):
        p = x[self.tets]
        return np.stack([p[:, i] - p[:, 0] for i in [1, 2, 3]], axis=2)

    def evaluate(self, x, bulk=30.0, gradient=True):
        F = self.edges(x) @ self.inverse
        J = np.linalg.det(F)
        if not np.all(np.isfinite(J)) or np.min(J) <= 1e-9:
            return np.inf, None, J
        Ic = np.einsum('tij,tij->t', F, F)
        logJ = np.log(J)
        scale = J ** (-2 / 3)
        energy = float(np.dot(self.weights, .5 * (scale * Ic - 3) + .5 * bulk * logJ**2))
        if not gradient:
            return energy, None, J
        inverseT = np.linalg.inv(F).transpose(0, 2, 1)
        P = scale[:, None, None] * (F - Ic[:, None, None] * inverseT / 3) + bulk * logJ[:, None, None] * inverseT
        G = (P @ self.inverse.transpose(0, 2, 1)) * self.weights[:, None, None]
        g = np.zeros_like(x)
        np.add.at(g, self.tets[:, 0], -G.sum(2))
        for k in range(3):
            np.add.at(g, self.tets[:, k + 1], G[:, :, k])
        return energy, g, J

    def predict_boundary(self, x, target):
        result = x.copy()
        delta = target - x[self.fixed]
        result[self.free] += self.factor.solve(-(self.coupling @ delta))
        result[self.fixed] = target
        return result

    def minimize(self, initial, bulk=30.0, iterations=200, residual_tolerance_mm=.001):
        """Preconditioned L-BFGS with a positive-J Armijo line search.
        Positions are in source millimetres. Quasi-static, no active force.
        An iteration limit is explicitly reported, never marked converged.
        """
        x = initial.copy()
        energy, fullg, J = self.evaluate(x, bulk)
        if not np.isfinite(energy):
            raise ValueError('Initial condition contains non-positive Jacobians')
        g = fullg[self.free]
        history = []
        gamma = 1 / (1 + bulk)
        reason = 'iteration-limit'
        accepted = 0
        for k in range(iterations):
            q = g.copy()
            alphas = []
            for s, y, rho in reversed(history):
                alpha = rho * np.sum(s * q)
                alphas.append(alpha)
                q -= alpha * y
            z = gamma * self.factor.solve(q)
            for (s, y, rho), alpha in zip(history, reversed(alphas)):
                z += s * (alpha - rho * np.sum(y * z))
            d = -z
            slope = float(np.sum(g * d))
            if slope >= 0:
                history = []
                d = -self.factor.solve(g) / (1 + bulk)
                slope = float(np.sum(g * d))
            max_direction = float(np.linalg.norm(d, axis=1).max(initial=0))
            if max_direction < residual_tolerance_mm:
                # A small L-BFGS scaling must not manufacture convergence.
                # Stopping uses the unscaled source-FEM metric; multiplying
                # the bulk penalty or changing gamma cannot lower this test.
                residual = self.factor.solve(g)
                if np.linalg.norm(residual, axis=1).max(initial=0) < residual_tolerance_mm:
                    reason = 'residual-tolerance'
                    break
                history = []
                gamma = 1 / (1 + bulk)
                d = -residual / (1 + bulk)
                slope = float(np.sum(g * d))
            trial = x.copy()
            gain = 1.0
            success = False
            for _ in range(40):
                trial[self.free] = x[self.free] + gain * d
                value, _, nextJ = self.evaluate(trial, bulk, gradient=False)
                if np.isfinite(value) and value <= energy + .0001 * gain * slope:
                    success = True
                    break
                gain *= .5
            if not success:
                reason = 'line-search-stalled'
                break
            nextE, nextg, J = self.evaluate(trial, bulk)
            nextg = nextg[self.free]
            s = gain * d
            y = nextg - g
            sy = float(np.sum(s * y))
            if sy > 1e-12 * np.linalg.norm(s) * np.linalg.norm(y):
                Hy = self.factor.solve(y)
                denominator = float(np.sum(y * Hy))
                if denominator > 0:
                    gamma = sy / denominator
                    history.append((s, y, 1 / sy))
                    history = history[-10:]
            x, energy, g = trial, nextE, nextg
            accepted += 1
        _, final_gradient, J = self.evaluate(x, bulk)
        residual = self.factor.solve(final_gradient[self.free])
        return x, {
            'reason': reason, 'converged': reason == 'residual-tolerance',
            'iterations': accepted, 'energy': energy, 'bulkPenalty': bulk,
            'minJ': float(J.min()), 'maxJ': float(J.max()),
            'JQuantiles': np.quantile(J, [0, .01, .5, .99, 1]).tolist(),
            'volumeRatio': float(np.dot(self.weights, J)),
            'nonPositiveElements': int((J <= 0).sum()),
            'preconditionedGradientMaxMm': float(np.linalg.norm(residual, axis=1).max(initial=0)),
            'convergenceMetric': 'source FEM inverse times gradient, without bulk or L-BFGS scaling',
            'residualToleranceMm': residual_tolerance_mm,
            'activation': None, 'force': None, 'boneContactEvaluated': False,
        }
