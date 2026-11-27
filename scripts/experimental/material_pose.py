"""Exact mappings between a source-mm material chart and final rigid bones.
No joint interpolation or anatomical angle model lives in this module.
"""
import numpy as np
from scipy.spatial.transform import Rotation


def rigid(frame):
    R = Rotation.from_quat(frame['q']).as_matrix()
    t = np.asarray(frame['p']) - R @ np.asarray(frame['rest'])
    return R, t


def material_to_bone(target, key, center, scale):
    """Row-vector x @ R + t, material mm -> rest bone chart mm."""
    RA, ta = rigid(target['a'])
    RB, tb = rigid(target[key])
    return RA.T @ RB, ((np.asarray(center) @ RA.T + ta - tb) @ RB - center) / scale


def material_to_world(nodes, target, center, scale):
    R, t = rigid(target['a'])
    return (np.asarray(nodes) * scale + center) @ R.T + t


def boundary_targets(rest, fixed, classes, target, center, scale, bone_keys=None):
    """Rigid attachment regions may refer to any number of final bone frames.
    Default two-bone data remain compatible; unregistered labels are rejected
    rather than silently following the material-reference bone.
    """
    result = np.asarray(rest)[fixed].copy()
    labels = np.asarray(classes)[fixed]
    keys = {0: 'a', 1: 'b'} if bone_keys is None else bone_keys
    unique = np.unique(labels)
    for label in unique:
        if label not in keys:
            raise ValueError('Attachment class has no registered bone: ' + str(label))
    for label in unique:
        key = keys[label]
        if key == 'a':
            continue
        R, t = material_to_bone(target, key, center, scale)
        mask = labels == label
        result[mask] = (result[mask] - t) @ R.T
    return result
