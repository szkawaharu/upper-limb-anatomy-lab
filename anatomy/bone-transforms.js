(function (root) {
  "use strict";
  // A snapshot can only be published after joint constraints have finished.
  // Bone-local coordinates use each transform's explicit rest pivot.
  class BoneSnapshot {
    constructor(transforms, revision) {
      if (!Number.isInteger(revision) || revision < 0)
        throw Error("Invalid skeleton revision");
      this.revision = revision;
      this.transforms = {};
      for (const [id, a] of Object.entries(transforms)) {
        if (
          a.q.length !== 4 ||
          a.p.length !== 3 ||
          a.rest.length !== 3 ||
          ![...a.q, ...a.p, ...a.rest].every(Number.isFinite) ||
          Math.abs(Math.hypot(...a.q) - 1) > 1e-5
        )
          throw Error("Invalid rigid transform " + id);
        this.transforms[id] = Object.freeze({
          q: Object.freeze([...a.q]),
          p: Object.freeze([...a.p]),
          rest: Object.freeze([...a.rest]),
        });
      }
      Object.freeze(this.transforms);
      Object.freeze(this);
    }
    get(id) {
      const a = this.transforms[id];
      if (!a) throw Error("Unknown bone " + id);
      return a;
    }
    point(id, local) {
      const a = this.get(id);
      return root.AnatomyMath.add(a.p, root.AnatomyMath.qrot(a.q, local));
    }
  }
  root.BoneSnapshot = BoneSnapshot;
})(globalThis);
