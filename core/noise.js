/*
 * core/noise.js — seeded PRNG, integer hashes, value noise and fBm in 2D and 3D.
 *
 * Pure functions, no DOM. Ported from legacy/terrain_sandbox/terrain-model.js
 * (mulberry32, hash2, valueNoise, fbm) and extended to 3D for lithology blobs
 * (core/lithology.js). `valueNoise` / `fbm` keep their legacy 2D names so d3
 * migrates as a drop-in.
 *
 * Exposed as GS.noise (browser) or module.exports (node). Mirrored in py/gs_core.py.
 */
(function (root) {
  "use strict";
  var GS = (root.GS = root.GS || {});

  // Seeded PRNG in [0, 1).
  function mulberry32(a) {
    return function () {
      a |= 0;
      a = (a + 0x6d2b79f5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // Deterministic [0, 1) at integer lattice point for a seed.
  function hash2(ix, iy, seed) {
    var h = Math.imul(ix, 374761393) + Math.imul(iy, 668265263) + Math.imul(seed, 982451653);
    h = Math.imul(h ^ (h >>> 13), 1274126177);
    h = h ^ (h >>> 16);
    return (h >>> 0) / 4294967296;
  }
  function hash3(ix, iy, iz, seed) {
    var h = Math.imul(ix, 374761393) + Math.imul(iy, 668265263) + Math.imul(iz, 1103515245) + Math.imul(seed, 982451653);
    h = Math.imul(h ^ (h >>> 13), 1274126177);
    h = h ^ (h >>> 16);
    return (h >>> 0) / 4294967296;
  }

  function smoothstep(t) { return t * t * (3 - 2 * t); }

  // Single-octave value noise at (x, y) in lattice units → [0, 1].
  function valueNoise2(x, y, seed) {
    var x0 = Math.floor(x), y0 = Math.floor(y);
    var fx = smoothstep(x - x0), fy = smoothstep(y - y0);
    var v00 = hash2(x0, y0, seed), v10 = hash2(x0 + 1, y0, seed);
    var v01 = hash2(x0, y0 + 1, seed), v11 = hash2(x0 + 1, y0 + 1, seed);
    var a = v00 + (v10 - v00) * fx;
    var b = v01 + (v11 - v01) * fx;
    return a + (b - a) * fy;
  }
  function valueNoise3(x, y, z, seed) {
    var x0 = Math.floor(x), y0 = Math.floor(y), z0 = Math.floor(z);
    var fx = smoothstep(x - x0), fy = smoothstep(y - y0), fz = smoothstep(z - z0);
    var v000 = hash3(x0, y0, z0, seed), v100 = hash3(x0 + 1, y0, z0, seed);
    var v010 = hash3(x0, y0 + 1, z0, seed), v110 = hash3(x0 + 1, y0 + 1, z0, seed);
    var v001 = hash3(x0, y0, z0 + 1, seed), v101 = hash3(x0 + 1, y0, z0 + 1, seed);
    var v011 = hash3(x0, y0 + 1, z0 + 1, seed), v111 = hash3(x0 + 1, y0 + 1, z0 + 1, seed);
    var a0 = v000 + (v100 - v000) * fx, b0 = v010 + (v110 - v010) * fx;
    var a1 = v001 + (v101 - v001) * fx, b1 = v011 + (v111 - v011) * fx;
    var c0 = a0 + (b0 - a0) * fy, c1 = a1 + (b1 - a1) * fy;
    return c0 + (c1 - c0) * fz;
  }

  // Fractal Brownian motion (summed octaves) → [0, 1].
  function fbm2(x, y, seed, octaves, persistence, lacunarity) {
    var amp = 1, freq = 1, sum = 0, norm = 0;
    for (var o = 0; o < octaves; o++) {
      sum += amp * valueNoise2(x * freq, y * freq, seed + o * 1013);
      norm += amp;
      amp *= persistence;
      freq *= lacunarity;
    }
    return sum / norm;
  }
  function fbm3(x, y, z, seed, octaves, persistence, lacunarity) {
    var amp = 1, freq = 1, sum = 0, norm = 0;
    for (var o = 0; o < octaves; o++) {
      sum += amp * valueNoise3(x * freq, y * freq, z * freq, seed + o * 1013);
      norm += amp;
      amp *= persistence;
      freq *= lacunarity;
    }
    return sum / norm;
  }

  var api = {
    mulberry32: mulberry32,
    hash2: hash2, hash3: hash3,
    smoothstep: smoothstep,
    valueNoise: valueNoise2, valueNoise2: valueNoise2, valueNoise3: valueNoise3,
    fbm: fbm2, fbm2: fbm2, fbm3: fbm3
  };
  GS.noise = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window !== "undefined" ? window : this);
