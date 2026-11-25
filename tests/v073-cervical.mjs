import assert from "node:assert/strict";
import "../engine/math.js";
import "../engine/cervical-kinematics.js";
const C = globalThis.CervicalKinematics,
  M = globalThis.AnatomyMath;
let count = 0;
for (let angle = -75; angle <= 75; angle++) {
  const p = C.pose(angle),
    s = C.segments;
  assert(
    Math.abs(Object.values(p.relative).reduce((a, b) => a + b, 0) - angle) <
      1e-9,
  );
  assert(Math.abs(p.relative.C1 / (angle || 1)) <= 0.60000001);
  for (const c of s) {
    const a = p.transforms[c.parent],
      b = p.transforms[c.id];
    assert(
      M.len(M.sub(M.apply(a, b.rest), b.p)) < 1e-9,
      "joint pivot continuity",
    );
    assert(Math.abs(Math.hypot(...b.q) - 1) < 1e-9);
  }
  assert(
    Math.abs(
      (2 * Math.atan2(p.transforms.C0.q[1], p.transforms.C0.q[3]) * 180) /
        Math.PI -
        angle,
    ) < 1e-8,
  );
  count++;
}
for (const yaw of [0, Math.PI / 2, Math.PI, Math.PI * 1.5])
  for (const angle of [-60, 0, 60]) {
    const d = C.beginDrag(20, 30, angle, yaw);
    assert.equal(C.drag(d, 20, 30), angle);
    assert.equal(C.drag(d, 20, 200), angle);
    assert.equal(C.drag(d, 50, 30), C.drag(d, 50, 30), "held input idempotent");
    assert(Object.isFrozen(d));
    assert(C.drag(d, 20000, 30) <= 75);
  }
const p = C.pose(0);
for (const x of Object.values(p.transforms))
  assert.deepEqual(x.q, [0, 0, 0, 1]);
assert.throws(() => C.pose(NaN));
console.log({
  poses: count,
  chain: C.segments.map((x) => x.id),
  C1C2Share: 0.6,
  range: C.range,
});
