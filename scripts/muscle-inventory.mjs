import { readFileSync, writeFileSync } from "node:fs";
const html = readFileSync("index.html", "utf8"),
  js = html.match(/<script>([\s\S]*?)<\/script>/)[1];
const source =
  js.slice(
    js.indexOf("const V43_CONFIG="),
    js.indexOf("const MEDICAL_ANATOMY="),
  ) + js.slice(js.indexOf("const V58_DEEP="), js.indexOf("const V60_DELT="));
const { muscles, tendons, atlas, cuff } = new Function(
  source +
    ";return {muscles:MUSCLE_PARTS,tendons:TENDON_PARTS,atlas:V60_ALL,cuff:SCAPULAR_MUSCLES}",
)();
const rows = atlas.map((s) => ({
  id: s.id,
  name: s.name,
  source: s.file,
  representation: "atlas material map",
  migration: "shared CPU/GPU kernel; attachments estimated",
  validation: "engineering sweep; clinical/contact unverified",
}));
for (const id of muscles)
  rows.push({
    id,
    source: "procedural",
    representation: cuff.includes(id)
      ? "curved sheet loft"
      : id === "pronatorQuadratus"
        ? "quadrilateral sheet"
        : "shaped multi-station path",
    migration:
      "adapter; legacy routing retained, not a complete material-map migration",
    validation:
      "legacy synthetic continuity tests; geometric contact unverified",
  });
for (const id of tendons)
  rows.push({
    id,
    source: "procedural",
    representation: "tendon/aponeurosis path",
    migration: "adapter; prescribed distal taper and routing retained",
    validation: "endpoint continuity only; sliding/contact unverified",
  });
writeFileSync(
  "docs/validation/v073/migration.json",
  JSON.stringify(rows, null, 2),
);
console.log({
  atlas: atlas.length,
  proceduralMuscles: muscles.length,
  proceduralTendons: tendons.length,
  total: rows.length,
});
