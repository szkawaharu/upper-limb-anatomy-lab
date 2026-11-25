# v0.73 implementation work log

Base: `025b768` (v0.72.0), preserved at `archive/v0.72.0-before-common-engine`.
Working branch: `feature/v0.73.0-common-engine`.

## Audit before edits

- Input → `compute` / `makeRealPose` → `shoulderPoseFor` → v0.66 thorax surface guide → v0.67 orientation/GH guard → v0.72 contact guard → final clavicle/scapula/humerus transforms → elbow/forearm/hand → procedural paths. Atlas muscle fields are generated later in `v58DrawDeep`, using these final transforms. Inactive arm uses a separate cached snapshot; any new pose dependency must invalidate that cache.
- 17 atlas muscle subdivisions use clustered BodyParts3D surfaces (1 mm; serratus/latissimus 1.5 mm), landmark/nearest-bone attachment weights, mesh adjacency smoothing, and multi-pose weight selection. Source frame is millimetres/Z-up; `atlasToWorld` is a uniform .00555 scale and axis permutation/reflection. Source commit `2c83c698788cff57f6317bddcf712e7acfbea7d8` is pinned for torso/muscles. Older arm sources still use main and need source pinning.
- Atlas path: LBS → optional DQS → v0.62 route/guard → bounded displacement → v0.64 fiber OR v0.65 9-slip OR v0.70 fan → v0.69 soft tissue (pectoralis only). v0.68.2 adjusts controls on non-fan muscles. CPU audit has corresponding point evaluators, but actual CPU/GPU readback comparison was absent. GPU failure previously falls back to an older method, rather than the same algorithm on CPU.
- Most arm/forearm/hand bellies and tendons are procedural Catmull-Rom/linear paths rebuilt on pose changes, not atlas surfaces. Rotator cuff/teres use distinct curved sheet lofts. Full lists are `MUSCLE_PARTS`, `TENDON_PARTS`, `V60_ALL`, `V43_CONFIG`. Individual heads are separate meshes. They must not be reported as measured fascicles or atlas muscles.
- Broad pectoral origin → humerus: v0.70 measures raw longitudinal/transverse derivatives but transports axes using blended bone rotation, not actual sheet tangent. Independent stretch clamps and a downstream rest-relative ellipsoid guard can leave inconsistent offsets/thickness under elevation. This is a geometric diagnosis, not a verified sole cause. New engine must preserve rest geometry, use an aligned stable material frame, and audit after contact.
- No activation input/force model exists. Existing elbow/wrist-angle-dependent visual bulge is an educational exaggeration, not measured activity. New engine must distinguish this from activation.
- STL weld/clustering can leave open/nonmanifold meshes. Only oriented closed manifold geometry supports signed volume as a material-volume proxy. Surface area is not a volumetric Jacobian. Baseline/new intersection must be distinguished; exact inter-muscle collision remains unverified until measured.
- Existing 21 offline test entry points were executed before changes. 15 pass; 4 pure-extraction harnesses fail due to outdated missing guide dependencies; 2 integration entry points exit with unsettled await because their mocked timers never resolve the newer loader's yield. See `validation/v073/baseline-tests.json`. These failures must be repaired without relaxing geometric thresholds.

## Execution plan / acceptance

1. Save baseline actual atlas poses/cameras/CPU rendering timing and source mesh topology.
2. Independent anatomy definitions + immutable bone snapshot + common path/frame/volume/contact/tendon/audit modules. Reuse regional fan/slip precomputation, not legacy postprocessing. Exclusive old/new dispatch; CPU same-algorithm fallback and recorded reasons.
3. Integrate atlas subdivisions and procedural adapters; maintain morphologies and endpoint identities. Pectoralis priority, then deltoid/scapular muscles. Document partial adapters accurately.
4. T1→C7→…→C1→C0 hierarchy, pinned matching skull parts and confirmed neck muscles; degree distribution is configurable educational approximation. Slider/reset and camera-locked direct drag with arm independence.
5. Repair baseline harness defects; run all legacy/new tests, browser actual-atlas sweeps, screenshots and GPU readback. Record max/outliers, null for unsupported measures, mobile device untested.
6. Commit each functional unit. Main/Pages release only if required scope is complete and no major regressions remain. Otherwise push work branch and write exact resume steps. No medical accuracy certification from engineering tests.
