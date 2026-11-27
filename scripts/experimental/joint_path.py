"""Transactional traversal of a frozen, integer-indexed actual joint-pose grid.
Rejected solves never change the accepted mesh or pose. Subdivision samples
actual precomputed joint frames, not SLERP of endpoint bones or boundary pins.
"""
import numpy as np


def follow(initial, start_index, destinations, solve_at, on_accept=None):
    x = np.array(initial, copy=True)
    index = int(start_index)
    accepted, rejected = [], []
    for requested in destinations:
        pending = [int(requested)]
        while pending:
            target = pending.pop(0)
            if target == index:
                continue
            trial, report = solve_at(x.copy(), index, target)
            if report.get('converged') and np.isfinite(trial).all():
                previous = index
                x, index = np.array(trial, copy=True), target
                record = {'fromIndex': previous, 'index': index, 'requestedIndex': int(requested), **report}
                accepted.append(record)
                if on_accept:
                    on_accept(x.copy(), index, accepted, rejected)
            else:
                rejected.append({'fromIndex': index, 'attemptedIndex': target, 'requestedIndex': int(requested), **report})
                distance = abs(target - index)
                if distance <= 1:
                    return x, {'completed': False, 'index': index, 'accepted': accepted, 'rejected': rejected}
                middle = index + (1 if target > index else -1) * (distance // 2)
                pending = [middle, target, *pending]
    return x, {'completed': True, 'index': index, 'accepted': accepted, 'rejected': rejected}
