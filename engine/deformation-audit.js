(function (root) {
  "use strict";
  const { sub, cross, dot, len } = root.AnatomyMath;
  function topology(positions, faces) {
    const edges = new Map();
    let degenerate = 0;
    for (let i = 0; i < faces.length; i += 3) {
      const ids = [faces[i], faces[i + 1], faces[i + 2]];
      if (ids.some((k) => k < 0 || k * 3 + 2 >= positions.length))
        throw Error("Invalid face");
      const p = ids.map((k) => Array.from(positions.slice(k * 3, k * 3 + 3)));
      if (len(cross(sub(p[1], p[0]), sub(p[2], p[0]))) < 1e-12) degenerate++;
      for (let j = 0; j < 3; j++) {
        const a = ids[j],
          b = ids[(j + 1) % 3],
          key = Math.min(a, b) + ":" + Math.max(a, b),
          e = edges.get(key) || [0, 0];
        e[0]++;
        e[1] += a < b ? 1 : -1;
        edges.set(key, e);
      }
    }
    let boundary = 0,
      nonmanifold = 0,
      inconsistent = 0;
    for (const [n, s] of edges.values()) {
      if (n === 1) boundary++;
      else if (n !== 2) nonmanifold++;
      if (n === 2 && s !== 0) inconsistent++;
    }
    return {
      closed: !boundary && !nonmanifold && !inconsistent && !degenerate,
      boundary,
      nonmanifold,
      inconsistent,
      degenerate,
    };
  }
  function signedVolume(p, f) {
    let v = 0;
    const o = p.slice(0, 3);
    for (let i = 0; i < f.length; i += 3) {
      const a = sub(p.slice(f[i] * 3, f[i] * 3 + 3), o),
        b = sub(p.slice(f[i + 1] * 3, f[i + 1] * 3 + 3), o),
        c = sub(p.slice(f[i + 2] * 3, f[i + 2] * 3 + 3), o);
      v += dot(a, cross(b, c)) / 6;
    }
    return v;
  }
  function audit(state, now, previous = null) {
    const B = state.binding,
      rest = B.rest,
      faces = B.faces,
      scale = B.definition.source.worldUnitsPerMm,
      t = B.topology || (B.topology = topology(rest, faces));
    let nonfinite = 0,
      attachmentError = 0,
      attachmentSamples = 0,
      maxStep = 0,
      minArea = Infinity,
      maxArea = 0,
      minEdge = Infinity,
      maxEdge = 0,
      degenerate = 0,
      minStretch = Infinity,
      maxStretch = 0,
      minDet = Infinity,
      maxDet = -Infinity,
      bounded = 0,
      invalidFrames = 0,
      conditionedFrames = 0,
      minEffectiveTangentDot = Infinity;
    for (let i = 0; i < B.weights.length; i++) {
      const p = now.slice(3 * i, 3 * i + 3),
        r = rest.slice(3 * i, 3 * i + 3),
        w = B.weights[i];
      if (![...p].every(Number.isFinite)) nonfinite++;
      if (w === 0 || w === 1) {
        attachmentSamples++;
        const a = w ? state.b : state.originAt ? state.originAt(r) : state.a;
        attachmentError = Math.max(
          attachmentError,
          len(sub(p, root.AnatomyMath.apply(a, r))) / scale,
        );
      }
      if (previous)
        maxStep = Math.max(
          maxStep,
          len(sub(p, previous.slice(3 * i, 3 * i + 3))) / scale,
        );
      const e = root.MuscleDeformation.evaluate(
        state,
        r,
        w,
        B.lane[i],
        null,
        B.material[i],
      );
      if (e.invalidFrame) invalidFrames++;
      if (e.frameConfidence != null && e.frameConfidence < 1) conditionedFrames++;
      if (e.effectiveTangentDot != null) minEffectiveTangentDot = Math.min(minEffectiveTangentDot,e.effectiveTangentDot);
      if (e.bounded) bounded++;
      if (e.stretch != null) {
        minStretch = Math.min(minStretch, e.stretch);
        maxStretch = Math.max(maxStretch, e.stretch);
        minDet = Math.min(minDet, e.determinant);
        maxDet = Math.max(maxDet, e.determinant);
      }
    }
    for (let i = 0; i < faces.length; i += 3) {
      const ids = [faces[i], faces[i + 1], faces[i + 2]],
        r = ids.map((k) => rest.slice(3 * k, 3 * k + 3)),
        p = ids.map((k) => now.slice(3 * k, 3 * k + 3)),
        a0 = len(cross(sub(r[1], r[0]), sub(r[2], r[0]))),
        a = len(cross(sub(p[1], p[0]), sub(p[2], p[0])));
      if (a < 1e-12) degenerate++;
      if (a0 > 1e-12) {
        minArea = Math.min(minArea, a / a0);
        maxArea = Math.max(maxArea, a / a0);
      }
      for (let j = 0; j < 3; j++) {
        const d0 = len(sub(r[j], r[(j + 1) % 3]));
        if (d0 > 1e-9) {
          const q = len(sub(p[j], p[(j + 1) % 3])) / d0;
          minEdge = Math.min(minEdge, q);
          maxEdge = Math.max(maxEdge, q);
        }
      }
    }
    const v0 = t.closed ? signedVolume(rest, faces) : null,
      v = t.closed ? signedVolume(now, faces) : null,
      validVolume = v0 !== null && Math.abs(v0) > 1e-12;
    return {
      id: B.definition.id,
      topology: t,
      nonfinite,
      attachmentMaxMm: attachmentSamples ? attachmentError : null,
      attachmentSamples,
      maxStepMm: previous ? maxStep : null,
      restVolume: validVolume ? Math.abs(v0) : null,
      volume: validVolume ? Math.abs(v) : null,
      volumeRatio: validVolume ? v / v0 : null,
      volumeStatus: validVolume ? "closed-oriented-proxy" : "not-evaluable",
      areaRatio: [minArea, maxArea],
      edgeRatio: [minEdge, maxEdge],
      newDegenerateFaces: Math.max(0, degenerate - t.degenerate),
      stretch: [minStretch, maxStretch],
      localDeterminantApprox: [minDet, maxDet],
      boundedVertices: bounded,
      invalidFrames,
      conditionedFrames,
      minEffectiveTangentDot:Number.isFinite(minEffectiveTangentDot)?minEffectiveTangentDot:null,
      volumetricJacobian: null,
      newIntersections: null,
      newInversions: null,
      boneDistance: null,
      clinicalValidation: false,
    };
  }
  root.DeformationAudit = Object.freeze({ topology, signedVolume, audit });
})(globalThis);
