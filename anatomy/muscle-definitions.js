/* Definition data is separate from deformation. Coordinates/architecture are
 * atlas-derived estimates unless explicitly marked otherwise. */
(function (root) {
  "use strict";
  const source = {
    atlas: "BodyParts3D",
    commit: "2c83c698788cff57f6317bddcf712e7acfbea7d8",
    license: "CC BY-SA 2.1 Japan",
    sourceUnits: "mm",
    worldUnitsPerMm: 0.00555,
  };
  const neckOrigins = Object.freeze({
    v58Levator: Object.freeze(["C1", "C2", "C3", "C4"]),
    v59TrapUpper: Object.freeze(["C0", "C2", "C4", "C7", "T1"]),
  });
  const configs = {
    v58RhomboidMajor: ["sheet", "thorax"],
    v58RhomboidMinor: ["sheet", "thorax"],
    v58Levator: ["fusiform", null],
    v58Serratus: ["slips", "thorax"],
    v58Subscap: ["sheet", "head"],
    v59TrapUpper: ["sheet", null],
    v59TrapMiddle: ["sheet", "thorax"],
    v59TrapLower: ["sheet", "thorax"],
    v59Lat: ["fan", "thorax"],
    v60PecClav: ["fan", "thorax"],
    v60PecStern: ["fan", "thorax"],
    v60PecAbd: ["fan", "thorax"],
    v60PecMinor: ["fan", "thorax"],
    v60Subclavius: ["fusiform", null],
    v60DeltAnt: ["pennate", "head"],
    v60DeltMid: ["pennate", "head"],
    v60DeltPost: ["pennate", "head"],
  };
  function create(spec, vertices, faces, weights, anchors) {
    const c = configs[spec.id] || [
        spec.morphology || "fusiform",
        spec.contact || null,
      ],
      fixed = (end) => {
        const bone = end ? spec.b : spec.a,
          a = anchors[end],
          indices = [],
          local = [];
        for (let i = 0; i < weights.length; i++)
          if (weights[i] === (end ? 1 : 0)) {
            indices.push(i);
            local.push(vertices[i].map((x, k) => x - a.rest[k]));
          }
        return {
          bone,
          kind: "surface-region",
          indices: Uint32Array.from(indices),
          boneLocal: Float32Array.from(local.flat()),
          confidence: "estimated",
        };
      };
    return {
      id: spec.id,
      name: spec.name,
      side: spec.side || "right",
      region: spec.region || "upper-body",
      morphology: c[0],
      source: {
        ...source,
        file: spec.file || null,
        geometry: spec.file ? "atlas-surface" : "procedural-estimate",
      },
      restMesh: { vertices, faces },
      attachments: [fixed(0), fixed(1)],
      originBones: spec.originBones || neckOrigins[spec.id] || [spec.a],
      insertionBones: [spec.b],
      fibers: {
        arrangement: c[0],
        directions: null,
        restLengths: null,
        measured: false,
      },
      restCrossSection: null,
      restVolume: null,
      tissues: {
        belly: true,
        tendonStart:
          spec.tendonStart ?? (spec.id.startsWith("v60Pec") ? 0.76 : 0.86),
        aponeurosis: null,
      },
      viaPoints: [],
      obstacles: c[1] ? [{ kind: c[1], method: "rest-relative-envelope" }] : [],
      fascia:
        c[1] === "thorax"
          ? { kind: "rest-relative-envelope", measured: false }
          : null,
      sliding: { tissues: [], direction: null },
      constraints: {
        stretch: [0.35, 2.8],
        width: [0.4, 2.5],
        thickness: [0.4, 1.6],
      },
      activation: null,
      externalLoad: null,
      forceModel: null,
      confidence: "estimated",
      unverified: [
        "attachment anatomy",
        "fascicle trajectories",
        "exact bone contact",
        "inter-muscle sliding",
        "material volume",
        "active force",
        "clinical accuracy",
      ],
    };
  }
  root.MuscleDefinitions = Object.freeze({
    source,
    configs,
    neckOrigins,
    create,
  });
})(globalThis);
