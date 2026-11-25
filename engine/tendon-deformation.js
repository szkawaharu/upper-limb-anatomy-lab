(function (root) {
  "use strict";
  const { smooth, len, sub } = root.AnatomyMath;
  function fraction(t, regions) {
    if (!regions || regions.tendonStart == null) return 0;
    return smooth(
      (t - regions.tendonStart) / Math.max(1 - regions.tendonStart, 1e-6),
    );
  }
  function pathLength(path) {
    let d = 0;
    for (let i = 1; i < path.length; i++) d += len(sub(path[i], path[i - 1]));
    return d;
  }
  // Keep prescribed endpoints/routing stations and distal taper; no muscle-like
  // radial expansion for tendon. Passive length is distinct from activation.
  function metrics(rest, now) {
    const L0 = pathLength(rest),
      L = pathLength(now);
    return {
      restLength: L0,
      length: L,
      stretch: L0 > 1e-9 ? L / L0 : null,
      activation: null,
      force: null,
    };
  }
  root.TendonDeformation = Object.freeze({ fraction, pathLength, metrics });
})(globalThis);
