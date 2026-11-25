import assert from "node:assert/strict";
import "../muscle-fiber-v064.js";
import "../pectoral-fan-v070.js";
for (const file of [
  "math",
  "fiber-path",
  "volume-constraint",
  "contact-constraint",
  "tendon-deformation",
  "material-coordinates",
  "muscle-deformation",
  "deformation-audit",
  "gpu-deformation",
])
  await import("../engine/" + file + ".js");
await import("../anatomy/muscle-definitions.js");
await import("../anatomy/bone-transforms.js");
const {
  AnatomyMath: M,
  MuscleDeformation: E,
  DeformationAudit: A,
  MuscleDefinitions: D,
  VolumeConstraint: V,
} = globalThis;
const vertices = [],
  weights = [],
  faces = [];
for (let j = 0; j <= 20; j++)
  for (let k = 0; k <= 10; k++) {
    const t = j / 20,
      u = (k / 10 - 0.5) * (1 - t * 0.85);
    vertices.push([t, u, 0.02 * Math.sin(Math.PI * t) * (1 - u * u)]);
    weights.push(t);
  }
for (let j = 0; j < 20; j++)
  for (let k = 0; k < 10; k++) {
    const i = j * 11 + k;
    faces.push(i, i + 11, i + 1, i + 1, i + 11, i + 12);
  }
const a = M.identity(),
  restB = { q: [0, 0, 0, 1], p: [1, 0, 0], rest: [1, 0, 0] },
  def = D.create(
    { id: "syntheticFan", a: "trunk", b: "humerus", morphology: "fan" },
    vertices,
    faces,
    weights,
    [a, restB],
  ),
  binding = E.bind(def, weights);
assert(binding.sheet);
let maxStep = 0,
  previous = null;
for (let degrees = 0; degrees <= 180; degrees++) {
  const b = {
      ...restB,
      q: M.qaxis([0, 0, 1], (degrees * Math.PI) / 180),
      p: [1, 0.002 * degrees, 0],
    },
    state = E.pose(binding, a, b),
    now = E.deform(state),
    report = A.audit(state, now, previous);
  assert.equal(report.nonfinite, 0);
  assert(report.attachmentMaxMm < 0.001);
  assert.equal(report.invalidFrames, 0);
  assert(report.localDeterminantApprox[0] > 0);
  if (degrees === 0)
    assert(
      Math.max(...now.map((x, i) => Math.abs(x - binding.rest[i]))) < 1e-6,
    );
  if (previous) maxStep = Math.max(maxStep, report.maxStepMm);
  previous = now;
}
assert(maxStep < 10, "one degree bounded transition");
const neutral = E.deform(E.pose(binding, a, restB));
assert(
  Math.max(...neutral.map((x, i) => Math.abs(x - binding.rest[i]))) < 1e-6,
  "no neutral drift",
);
for (const shape of [true, false])
  for (const x of [0.4, 0.7, 1, 1.4, 2]) {
    const s = V.scales(x, 1.2, shape);
    if (!s.bounded) assert(Math.abs(s.determinant - 1) < 1e-10);
  }
const tet = new Float32Array([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1]),
  f = [0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];
assert(A.topology(tet, f).closed);
assert(!A.topology(tet, f.slice(0, 9)).closed);
assert(Math.abs(A.signedVolume(tet, f) - 1 / 6) < 1e-12);
const env = { center: [0, 0, 0], radii: [1, 1, 1], scale: 0.001 };
assert.deepEqual(
  globalThis.ContactConstraint.project([1, 0, 0], [2, 0, 0], 0, env),
  [2, 0, 0],
);
assert(
  M.len(globalThis.ContactConstraint.project([1, 0, 0], [2, 0, 0], 0.5, env)) <
    2,
);
assert.equal(globalThis.TendonDeformation.fraction(1, { tendonStart: 0.8 }), 1);
const tendon = V.scales(0.8, 1.2, true, 1);
assert.equal(tendon.u, 1);
assert.equal(tendon.n, 1);
assert.throws(
  () =>
    new globalThis.BoneSnapshot(
      { a: { q: [0, 0, 0, 2], p: [0, 0, 0], rest: [0, 0, 0] } },
      0,
    ),
);
assert.throws(() => new globalThis.BoneSnapshot({}, 0).get("missing"));
console.log(
  JSON.stringify({
    poses: 181,
    maxStepMm: maxStep,
    attachments: "< .001 mm",
    openMeshVolume: "not evaluable",
    neutralDrift: 0,
  }),
);

assert.throws(
  () =>
    E.bind(
      { ...def, restMesh: { vertices, faces: [0, 1, vertices.length] } },
      weights,
    ),
  /Invalid muscle faces/,
);
