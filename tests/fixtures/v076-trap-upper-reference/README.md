# Full-source numerical fixture — not anatomical acceptance

`solid.npz` is derived from BodyParts3D / Anatomography, © The Database Center for Life Science, under [CC BY-SA 2.1 Japan](https://creativecommons.org/licenses/by-sa/2.1/jp/). The source commit, original file, citation, modifications and repair audit are in `metadata.json`. The modified anatomical data retain that attribution and licence.

The original closed mesh has confirmed self-intersections. This fixture uses a separately audited repair candidate; components were not removed or joined. Source-to-candidate sampled vertex/surface distance is at most 0.1091mm; this is not a continuous Hausdorff certificate. Inferred attachment masks are not measured entheses.

The real source regression exercises the common offline solver at fixed 3°/6° abduction, with neck zero. It preserves repaired boundary geometry, exact bone attachment transforms and positive element Jacobians. The motion regression guard equals prescribed fixed-boundary displacement plus the measured source median edge. Global ±5% volume is the user's engineering target. Local J outliers remain and are reported. Passing does not qualify anatomy, the production app, full ROM, GPU solving or mobile performance.

Generation uses separate offline GPL MeshFix and AGPL TetGen tools. Neither library is a web runtime or CI dependency; CI consumes this immutable data fixture using numpy/scipy. Regeneration is documented in `docs/CONTACT_REFERENCE_V076.md`.
