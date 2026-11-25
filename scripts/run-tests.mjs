import { readdirSync, mkdirSync, writeFileSync } from "node:fs";
import { spawnSync } from "node:child_process";
const result = [];
mkdirSync("test-results", { recursive: true });
for (const file of readdirSync("tests")
  .filter(
    (x) =>
      x.endsWith(".mjs") && !x.startsWith("browser-") && x !== "v064-atlas.mjs",
  )
  .sort()) {
  const t = performance.now(),
    r = spawnSync(process.execPath, ["tests/" + file], {
      encoding: "utf8",
      timeout: 120000,
    }),
    output = (r.stdout || "") + (r.stderr || "");
  writeFileSync("test-results/" + file + ".log", output);
  result.push({
    test: file,
    code: r.status,
    error: r.error?.message,
    ms: Math.round(performance.now() - t),
  });
  console.log((r.status === 0 ? "PASS " : "FAIL ") + file);
  if (r.status !== 0) console.log(output.slice(-1800));
}
writeFileSync("test-results/results.json", JSON.stringify(result, null, 2));
process.exitCode = result.every((x) => x.code === 0) ? 0 : 1;
