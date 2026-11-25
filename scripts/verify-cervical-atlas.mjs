import { readFileSync, writeFileSync, mkdirSync } from "node:fs";
import { createHash } from "node:crypto";
import { join } from "node:path";
import "../anatomy/cervical-data.js";
const dir = process.argv[2];
if (!dir) throw Error("Supply a directory with the 28 pinned STL files");
const rows = [];
for (const e of globalThis.CervicalAtlasData) {
  const b = readFileSync(join(dir, e.file)),
    hash = createHash("sha1")
      .update(`blob ${b.length}\0`)
      .update(b)
      .digest("hex");
  if (b.length !== e.bytes || hash !== e.sha)
    throw Error("Source mismatch " + e.file);
  const n = b.readUInt32LE(80);
  if (b.length !== 84 + n * 50) throw Error("Invalid STL " + e.file);
  const lo = [Infinity, Infinity, Infinity],
    hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < n; i++)
    for (let j = 0; j < 3; j++)
      for (let k = 0; k < 3; k++) {
        const x = b.readFloatLE(96 + 50 * i + 12 * j + 4 * k);
        if (!Number.isFinite(x)) throw Error("Invalid position");
        lo[k] = Math.min(lo[k], x);
        hi[k] = Math.max(hi[k], x);
      }
  rows.push({
    id: e.id,
    name: e.name,
    kind: e.kind,
    sha: hash,
    bytes: b.length,
    triangles: n,
    atlasBoundsMm: [lo, hi],
  });
}
mkdirSync("test-results", { recursive: true });
writeFileSync(
  "test-results/cervical-source.json",
  JSON.stringify(rows, null, 2),
);
console.log(
  `Verified ${rows.length} STL byte counts and Git blob hashes at pinned source commit`,
);
