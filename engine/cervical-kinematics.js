(function (root) {
  "use strict";
  const M = root.AnatomyMath,
    { identity, apply, qaxis, qmul, clamp } = M;
  // Configuration is an educational interpolation, not an in-vivo measurement.
  // C1 denotes motion of atlas relative to axis (C1-C2).
  const segments = Object.freeze(
    [
      { id: "C7", parent: "T1", share: 0.02, pivot: [0, -72, 1391] },
      { id: "C6", parent: "C7", share: 0.04, pivot: [0, -77, 1404] },
      { id: "C5", parent: "C6", share: 0.07, pivot: [0, -81, 1416] },
      { id: "C4", parent: "C5", share: 0.08, pivot: [0, -84, 1430] },
      { id: "C3", parent: "C4", share: 0.08, pivot: [0, -85, 1442] },
      { id: "C2", parent: "C3", share: 0.07, pivot: [0, -85, 1455] },
      { id: "C1", parent: "C2", share: 0.6, pivot: [0, -83, 1480] },
      { id: "C0", parent: "C1", share: 0.04, pivot: [0, -81, 1491] },
    ].map((x) =>
      Object.freeze({
        ...x,
        axis: Object.freeze([0, 1, 0]),
        pivot: Object.freeze(x.pivot),
        confidence: "estimated-atlas-pivot",
      }),
    ),
  );
  function pose(degrees, toWorld = (x) => x, config = segments) {
    if (!Number.isFinite(degrees)) throw Error("Invalid neck angle");
    if (Math.abs(config.reduce((s, c) => s + c.share, 0) - 1) > 1e-9)
      throw Error("Cervical shares must sum to one");
    const angle = clamp(degrees, -75, 75),
      transforms = { T1: identity() },
      relative = {};
    for (const c of config) {
      const parent = transforms[c.parent];
      if (!parent || c.share < 0) throw Error("Invalid cervical hierarchy");
      const rest = toWorld(c.pivot),
        q = qmul(parent.q, qaxis(c.axis, (angle * c.share * Math.PI) / 180));
      transforms[c.id] = { q, p: apply(parent, rest), rest };
      relative[c.id] = angle * c.share;
    }
    return {
      angle,
      transforms,
      relative,
      flexion: 0,
      lateralBending: 0,
      clinicalValidation: false,
    };
  }
  function beginDrag(x, y, angle, yaw, pixelsPerDegree = 3) {
    const a = (angle * Math.PI) / 180,
      derivative = -Math.cos(yaw - a),
      direction =
        Math.abs(derivative) > 0.25
          ? Math.sign(derivative)
          : Math.sign(Math.sin(yaw - a)) || 1;
    return Object.freeze({ x, y, angle, yaw, direction, pixelsPerDegree });
  }
  function drag(start, x, y) {
    return clamp(
      start.angle + ((x - start.x) * start.direction) / start.pixelsPerDegree,
      -75,
      75,
    );
  }
  root.CervicalKinematics = Object.freeze({
    version: "0.73.0",
    segments,
    pose,
    beginDrag,
    drag,
    range: [-75, 75],
    distributionMeasured: false,
  });
})(globalThis);
