(function (root) {
  "use strict";
  const M = root.AnatomyMath,
    {
      add,
      sub,
      mul,
      dot,
      cross,
      len,
      norm,
      mix,
      qrot,
      qmix,
      apply,
      align,
      clamp,
      smooth,
    } = M;
  function curve(c, t) {
    const u = 1 - t;
    return add(
      add(mul(c[0], u * u * u), mul(c[1], 3 * u * u * t)),
      add(mul(c[2], 3 * u * t * t), mul(c[3], t * t * t)),
    );
  }
  function tangent(c, t) {
    const u = 1 - t;
    return add(
      add(mul(sub(c[1], c[0]), 3 * u * u), mul(sub(c[2], c[1]), 6 * u * t)),
      mul(sub(c[3], c[2]), 3 * t * t),
    );
  }
  function laneCurve(c, l) {
    const x = clamp(l, 0, c.length - 1),
      i = Math.min(c.length - 2, Math.floor(x));
    return c[i].map((p, j) => mix(p, c[i + 1][j], x - i));
  }
  function frame(c, l, t, reference) {
    const r = laneCurve(c, l),
      T = tangent(r, t),
      tl = len(T);
    if (tl < 1e-7) return null;
    const u = norm(T),
      l0 = clamp(l - 0.25, 0, c.length - 1),
      l1 = clamp(l + 0.25, 0, c.length - 1),
      U = mul(
        sub(curve(laneCurve(c, l1), t), curve(laneCurve(c, l0), t)),
        1 / Math.max(l1 - l0, 1e-6),
      );
    let v = sub(U, mul(u, dot(U, u)));
    if (len(v) < 1e-7) v = sub(reference, mul(u, dot(reference, u)));
    if (len(v) < 1e-7) return null;
    v = norm(v);
    return {
      c: curve(r, t),
      t: u,
      u: v,
      n: norm(cross(u, v)),
      tl,
      ul: len(sub(U, mul(u, dot(U, u)))),
    };
  }
  function pose(field, a, b) {
    return field.controls.map((c) =>
      c.map((p, j) => mix(apply(a, p), apply(b, p), smooth(j / 3))),
    );
  }
  function transported(R, S, a, b, t) {
    const q = qmix(a.q, b.q, smooth(t)),
      rt = qrot(q, R.t),
      ru = qrot(q, R.u),
      rn = qrot(q, R.n),
      tr = M.conditionedTransport(rt, S.t),
      effectiveT = qrot(tr.q, rt),
      u = qrot(tr.q, ru),
      n = qrot(tr.q, rn);
    return { ...S, rawT:S.t, t:effectiveT, u, n, frameConfidence:tr.confidence, inputTangentDot:tr.inputDot, effectiveTangentDot:dot(effectiveT,S.t) };
  }
  root.FiberPath = Object.freeze({
    curve,
    tangent,
    laneCurve,
    frame,
    pose,
    transported,
  });
})(globalThis);
