(function (root) {
  "use strict";
  const { clamp } = root.AnatomyMath;
  // Bounded local geometric approximation; hitting a bound is reported, never
  // called exact incompressibility. Sheet area uses perpendicular transverse width.
  function scales(longitudinal, transverse, sheet, tendon = 0) {
    const rawT = longitudinal,
      rawU = transverse,
      t = clamp(rawT, 0.35, 2.8),
      u = sheet ? clamp(rawU, 0.4, 2.5) : 1 / Math.sqrt(t);
    const n0 = sheet ? 1 / (t * u) : u,
      n = clamp(n0, 0.4, 1.6),
      k = clamp(tendon, 0, 1);
    const U = u * (1 - k) + k,
      N = n * (1 - k) + k;
    return {
      t,
      u: U,
      n: N,
      determinant: t * U * N,
      bounded: t !== rawT || (sheet && u !== rawU) || n !== n0,
    };
  }
  root.VolumeConstraint = Object.freeze({ scales });
})(globalThis);
