(function (root) {
  "use strict";
  class MuscleRenderer {
    constructor(gl, fragment, createMesh) {
      this.gl = gl;
      this.createMesh = createMesh;
      this.resources = new WeakMap();
      this.reason = null;
      this.backend = "gpu";
      try {
        const s = (kind, src) => {
          const x = gl.createShader(kind);
          gl.shaderSource(x, src);
          gl.compileShader(x);
          if (!gl.getShaderParameter(x, gl.COMPILE_STATUS))
            throw Error(gl.getShaderInfoLog(x));
          return x;
        };
        const vs = s(gl.VERTEX_SHADER, root.MuscleGPU.vertexShader()),
          fs = s(gl.FRAGMENT_SHADER, fragment),
          p = gl.createProgram();
        gl.attachShader(p, vs);
        gl.attachShader(p, fs);
        gl.linkProgram(p);
        gl.deleteShader(vs);
        gl.deleteShader(fs);
        if (!gl.getProgramParameter(p, gl.LINK_STATUS))
          throw Error(gl.getProgramInfoLog(p));
        this.program = p;
        this.attributes = Object.fromEntries(
          ["pos", "normal", "weight", "lane", "materialT", "originBlend"].map(
            (n) => [n, gl.getAttribLocation(p, n)],
          ),
        );
        this.uniforms = Object.fromEntries(
          [
            "mvp",
            "model",
            "anchorA",
            "anchorB",
            "originAnchor",
            "originQA",
            "originMixOn",
            "qa",
            "qb",
            "restControls[0]",
            "nowControls[0]",
            "lanes",
            "sheet",
            "tendonStart",
            "reference",
            "contactOn",
            "contactProximal",
            "unitScale",
            "contactIn",
            "contactOut",
            "contactMax",
            "contactCenter",
            "contactRestCenter",
            "contactRadii",
            "color",
            "opacity",
            "kind",
          ].map((n) => [n, gl.getUniformLocation(p, n)]),
        );
      } catch (e) {
        this.backend = "cpu";
        this.reason = String(e);
      }
    }
    resource(B) {
      let r = this.resources.get(B);
      if (r) return r;
      const gl = this.gl,
        buffer = (a) => {
          const b = gl.createBuffer();
          gl.bindBuffer(gl.ARRAY_BUFFER, b);
          gl.bufferData(gl.ARRAY_BUFFER, a, gl.STATIC_DRAW);
          return b;
        };
      r = {
        weight: buffer(B.weights),
        lane: buffer(B.lane),
        materialT: buffer(B.material),
        originBlend: buffer(
          B.originBlend || new Float32Array(B.weights.length),
        ),
        rest: new Float32Array(108),
        now: new Float32Array(108),
        cpu: new Map(),
      };
      this.resources.set(B, r);
      return r;
    }
    draw(state, g, matrices, color) {
      const gl = this.gl,
        B = state.binding,
        r = this.resource(B),
        u = this.uniforms;
      if (this.backend !== "gpu" || (state.originAt && !state.originMix))
        return false;
      gl.useProgram(this.program);
      r.rest.fill(0);
      r.now.fill(0);
      r.rest.set(state.restFlat);
      r.now.set(state.nowFlat);
      for (const [n, b, size] of [
        ["pos", g.vbo, 3],
        ["normal", g.nbo, 3],
        ["weight", r.weight, 1],
        ["lane", r.lane, 1],
        ["materialT", r.materialT, 1],
        ["originBlend", r.originBlend, 1],
      ]) {
        const a = this.attributes[n];
        if (a >= 0) {
          gl.bindBuffer(gl.ARRAY_BUFFER, b);
          gl.enableVertexAttribArray(a);
          gl.vertexAttribPointer(a, size, gl.FLOAT, false, 0, 0);
        }
      }
      for (const n of ["mvp", "model", "anchorA", "anchorB"])
        gl.uniformMatrix4fv(u[n], false, matrices[n]);
      gl.uniform1f(u.originMixOn, state.originMix ? 1 : 0);
      gl.uniform4fv(u.originQA, state.originMix?.q || state.a.q);
      gl.uniformMatrix4fv(
        u.originAnchor,
        false,
        matrices.originAnchor || matrices.anchorA,
      );
      gl.uniform4fv(u.qa, state.a.q);
      gl.uniform4fv(u.qb, state.b.q);
      gl.uniform3fv(u["restControls[0]"], r.rest);
      gl.uniform3fv(u["nowControls[0]"], r.now);
      gl.uniform1f(u.lanes, B.field.controls.length);
      gl.uniform1f(u.sheet, B.sheet ? 1 : 0);
      gl.uniform1f(u.tendonStart, B.definition.tissues.tendonStart);
      gl.uniform3fv(u.reference, B.reference);
      const c = state.contact;
      gl.uniform1f(u.contactOn, c ? 1 : 0);
      if (c) {
        gl.uniform1f(u.contactProximal, c.proximal ? 1 : 0);
        gl.uniform1f(u.unitScale, c.scale || 0.00555);
        gl.uniform1f(u.contactIn, c.inMm ?? 2);
        gl.uniform1f(u.contactOut, c.outMm ?? 6);
        gl.uniform1f(u.contactMax, c.maxMm ?? 20);
        gl.uniform3fv(u.contactCenter, c.center);
        gl.uniform3fv(u.contactRestCenter, c.restCenter || c.center);
        gl.uniform3fv(u.contactRadii, c.radii);
      }
      gl.uniform3fv(u.color, color);
      gl.uniform1f(u.opacity, 1);
      gl.uniform1f(u.kind, 1);
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, g.ibo);
      gl.drawElements(
        gl.TRIANGLES,
        g.count,
        g.indexType || gl.UNSIGNED_SHORT,
        0,
      );
      return true;
    }
    cpuGeometry(state, key) {
      const B = state.binding,
        r = this.resource(B);
      if (r.cpu.has(key)) return r.cpu.get(key);
      const p = root.MuscleDeformation.deform(state),
        vertices = Array.from({ length: p.length / 3 }, (_, i) =>
          Array.from(p.slice(i * 3, i * 3 + 3)),
        ),
        g = this.createMesh(vertices, Array.from(B.faces));
      if (r.cpu.size >= 2) {
        const old = r.cpu.values().next().value;
        for (const b of [old.vbo, old.nbo, old.ibo]) this.gl.deleteBuffer(b);
        r.cpu.delete(r.cpu.keys().next().value);
      }
      r.cpu.set(key, g);
      return g;
    }
  }
  root.MuscleRenderer = MuscleRenderer;
})(globalThis);
