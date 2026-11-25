// Curated static artifact: never publish the Git checkout, dependencies or logs.
import {
  cpSync,
  existsSync,
  mkdirSync,
  readFileSync,
  readdirSync,
  rmSync,
  writeFileSync,
} from "node:fs";
import { execFileSync } from "node:child_process";
import { resolve } from "node:path";

const output = resolve("_site");
const pkg = JSON.parse(readFileSync("package.json", "utf8"));
const commit = process.env.GITHUB_SHA || execFileSync("git", ["rev-parse", "HEAD"], { encoding: "utf8" }).trim();
if (!/^[0-9a-f]{40}$/.test(commit)) throw Error("Invalid deployment commit");

rmSync(output, { recursive: true, force: true });
mkdirSync(output, { recursive: true });
const publicFiles = readdirSync(".").filter((name) => /\.(html|js|zip)$/.test(name));
for (const name of [...publicFiles, "anatomy", "engine", "docs", "README.md", "HANDOFF_V073.md"]) {
  if (!existsSync(name)) throw Error("Missing public source: " + name);
  cpSync(name, resolve(output, name), { recursive: true });
}
writeFileSync(resolve(output, ".nojekyll"), "");
writeFileSync(resolve(output, "build-info.json"), JSON.stringify({
  version: pkg.version,
  commit,
  builtAt: new Date().toISOString(),
  validationStatus: "in-progress-with-known-issues",
  handoff: "HANDOFF_V073.md",
  knownIssues: "docs/KNOWN_ISSUES_V073.md",
}, null, 2) + "\n");

const html = readFileSync(resolve(output, "index.html"), "utf8");
if (!html.includes(`name="app-version" content="${pkg.version}"`)) throw Error("App/package version mismatch");
for (const match of html.matchAll(/<script\b[^>]*src="(\.[^"]+)"/g)) {
  const relative = match[1].split("?")[0];
  if (!existsSync(resolve(output, relative))) throw Error("Missing browser module: " + relative);
}
for (const excluded of [".git", ".github", "node_modules", "test-results"]) {
  if (existsSync(resolve(output, excluded))) throw Error("Private build input staged: " + excluded);
}
console.log(JSON.stringify({ version: pkg.version, commit, artifact: "_site", localScripts: [...html.matchAll(/<script\b[^>]*src="/g)].length }));
