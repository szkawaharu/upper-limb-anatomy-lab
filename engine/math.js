/* Shared geometric primitives; world units are declared by each anatomy source. */
(function (root) {
  "use strict";
  const clamp = (x, a, b) => Math.max(a, Math.min(b, x)),
    smooth = (x) => {
      x = clamp(x, 0, 1);
      return x * x * (3 - 2 * x);
    },
    add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]],
    sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]],
    mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s],
    dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2],
    cross = (a, b) => [
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0],
    ],
    len = (a) => Math.hypot(...a),
    norm = (a) => mul(a, 1 / Math.max(len(a), 1e-12)),
    mix = (a, b, t) => add(mul(a, 1 - t), mul(b, t));
  const qrot = (q, v) => {
      const t = mul(cross(q, v), 2);
      return add(v, add(mul(t, q[3]), cross(q, t)));
    },
    qmul = (a, b) => [
      a[3] * b[0] + a[0] * b[3] + a[1] * b[2] - a[2] * b[1],
      a[3] * b[1] - a[0] * b[2] + a[1] * b[3] + a[2] * b[0],
      a[3] * b[2] + a[0] * b[1] - a[1] * b[0] + a[2] * b[3],
      a[3] * b[3] - dot(a, b),
    ],
    qaxis = (axis, r) => [...mul(norm(axis), Math.sin(r / 2)), Math.cos(r / 2)],
    qmix = (a, b, t) => {
      const s = dot(a, b) + a[3] * b[3] < 0 ? -1 : 1,
        q = a.map((x, i) => x * (1 - t) + b[i] * s * t),
        n = Math.hypot(...q);
      return q.map((x) => x / Math.max(n, 1e-12));
    };
  const apply = (a, p) => add(a.p, qrot(a.q, sub(p, a.rest))),
    identity = () => ({ q: [0, 0, 0, 1], p: [0, 0, 0], rest: [0, 0, 0] });
  function align(u, v, p, axis) {
    const d = clamp(dot(u, v), -1, 1);
    if (d < -0.9999) return qrot(qaxis(axis, Math.PI), p);
    const k = cross(u, v),
      a = cross(k, p);
    return add(p, add(a, mul(cross(k, a), 1 / Math.max(1 + d, 1e-6))));
  }
  function conditionedTransport(u, v) {
    const d = clamp(dot(u, v), -1, 1), k = cross(u, v), s = len(k);
    // Engineering conditioning interval: Rodrigues' 1/(1+d)
    // grows from 100 to 10000 (the old hard-switch condition).
    const confidence = smooth(((1 + d) - 1e-4) / (1e-2 - 1e-4));
    const angle = confidence * Math.atan2(s, d);
    const q = s > 1e-12 && confidence > 0 ? qaxis(mul(k, 1 / s), angle) : [0,0,0,1];
    return {q, confidence, inputDot:d};
  }
  root.AnatomyMath = Object.freeze({
    clamp,
    smooth,
    add,
    sub,
    mul,
    dot,
    cross,
    len,
    norm,
    mix,
    qrot,
    qmul,
    qaxis,
    qmix,
    apply,
    identity,
    align,
    conditionedTransport,
  });
})(globalThis);
