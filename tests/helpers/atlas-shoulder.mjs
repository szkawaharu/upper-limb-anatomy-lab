// Reproduce the published app's shoulder inputs and thorax profile from the
// captured, pinned atlas fixture. No STL download or rendering in this helper.
import { readFileSync } from "node:fs";
import "../../engine/math.js";
import "../../scapulothoracic-v066.js";
import "../../scapulohumeral-v067.js";
import "../../scapulothoracic-v072.js";
export const fixture = JSON.parse(readFileSync(new URL("../fixtures/v073-continuity.json", import.meta.url)));
const M = globalThis.AnatomyMath;
const { add, sub, mul, qaxis, qmul, qrot, norm, clamp } = M;
const rad = x => x * Math.PI / 180;
const { anchors: A, data, scale } = fixture.context;
const cv = sub(A.AC, A.SC);
export function profileAt(p) {
  const { rows, dy } = data.profile;
  if (p[1] < rows[0].y - dy || p[1] > rows.at(-1).y + dy) return null;
  let a = rows[0], b = rows[0];
  for (let i = 1; i < rows.length; i++) { b = rows[i]; if (p[1] <= b.y) break; a = b; }
  if (p[1] > rows.at(-1).y) a = b = rows.at(-1);
  const t = a === b ? 0 : clamp((p[1] - a.y) / (b.y - a.y), 0, 1);
  return { c: a.c.map((x, i) => x + (b.c[i] - x) * t), r: a.r.map((x, i) => x + (b.r[i] - x) * t) };
}
export function gapFn(p) {
  const s = profileAt(p); if (!s) return null;
  const x = p[0] - s.c[0], z = p[2] - s.c[1], r = Math.hypot(x / s.r[0], z / s.r[1]);
  return r > 1e-8 ? Math.hypot(x, z) * (1 - 1 / r) : -Math.min(...s.r);
}
export function normalFn(p) {
  const s = profileAt(p); if (!s) return null;
  return norm([(p[0] - s.c[0]) / (s.r[0] ** 2), 0, (p[2] - s.c[1]) / (s.r[1] ** 2)]);
}
export function inputs(elev, plane, axial = 0) {
  const k = elev / 180, forward = Math.sin(rad(plane));
  const qC = qmul(qmul(qaxis([0,0,1], rad(14*k)), qaxis([0,1,0], rad(-(18-10*forward)*k))), qaxis(cv, rad(25*k)));
  const qS = qmul(qaxis([0,1,0], rad((-8+12*forward)*k)), qmul(qaxis([0,0,1], rad(elev/3)), qaxis([1,0,0], rad(13*k))));
  const qH = qmul(qaxis([Math.sin(rad(plane)),0,Math.cos(rad(plane))], rad(elev)), qaxis([0,-1,0], rad(axial)));
  const AC = add(A.SC, qrot(qC, cv));
  const restShaft=norm(sub(A.EL,A.GH)), guard=globalThis.ScapulohumeralGuardV067;
  const constraints=q=>guard.inequalities(guard.audit(q,qH,restShaft,elev,plane,qS));
  return {baseQC:qC, baseQS:qS, base:{qC,qS,AC}, SC:A.SC, cv, restAC:A.AC, elev, plane, qH, restShaft, points:data.v, gapFn, normalFn, scale, constraints};
}
export function solve(elev, plane, axial = 0) {
  const v66 = globalThis.ScapulothoracicEngineV066, v67 = globalThis.ScapulohumeralGuardV067, v72 = globalThis.ScapulothoracicContactGuardV072;
  const opts = inputs(elev,plane,axial), candidate = v66.solve(opts);
  const guarded = v67.refine({...opts, v66, candidate});
  const final = v72.refine({...opts, v66, v67, input:guarded});
  return {opts,candidate,guarded,final};
}
