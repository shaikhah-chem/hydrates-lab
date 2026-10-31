// Builds ball-and-stick GLB models of hydrated / anhydrous salts in the browser (no model files needed).
(function () {
  const COL = {Co: 'F090A0', Cl: '1FD01F', O: 'FF2020', H: 'F4F4F4', Cu: 'C88033', S: 'FFD830', Ba: '8E5CD9', bond: 'B0B8C0'};
  const RAD = {Co: 0.045, Cl: 0.05, O: 0.032, H: 0.019, Cu: 0.045, S: 0.04, Ba: 0.06};
  const lin = h => [0, 2, 4].map(i => { const x = parseInt(h.substr(i, 2), 16) / 255; return x <= 0.04045 ? x / 12.92 : Math.pow((x + 0.055) / 1.055, 2.4); });
  const add = (a, b) => a.map((v, i) => v + b[i]), sc = (a, s) => a.map(v => v * s);
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const norm = a => { const l = Math.hypot(...a); return a.map(v => v / l); };

  function sphere(c, r, lat = 16, lon = 24) {
    const v = [], n = [], idx = [];
    for (let i = 0; i <= lat; i++) { const th = Math.PI * i / lat;
      for (let j = 0; j <= lon; j++) { const ph = 2 * Math.PI * j / lon;
        const d = [Math.sin(th) * Math.cos(ph), Math.cos(th), Math.sin(th) * Math.sin(ph)];
        n.push(d); v.push(add(c, sc(d, r))); } }
    for (let i = 0; i < lat; i++) for (let j = 0; j < lon; j++) {
      const a = i * (lon + 1) + j, b = a + lon + 1; idx.push(a, a + 1, b, b, a + 1, b + 1); }
    return [v, n, idx];
  }
  function cylinder(p1, p2, r, seg = 14) {
    const ax = norm(p2.map((x, i) => x - p1[i]));
    const t = Math.abs(ax[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
    const u = norm(cross(ax, t)), w = cross(ax, u);
    const v = [], n = [], idx = [];
    for (let j = 0; j <= seg; j++) { const a = 2 * Math.PI * j / seg;
      const d = add(sc(u, Math.cos(a)), sc(w, Math.sin(a)));
      for (const p of [p1, p2]) { v.push(add(p, sc(d, r))); n.push(d); } }
    for (let j = 0; j < seg; j++) { const a = 2 * j; idx.push(a, a + 2, a + 1, a + 1, a + 2, a + 3); }
    return [v, n, idx];
  }
  function disc(y, r, h = 0.012, seg = 64) {
    const v = [[0, y + h, 0]], n = [[0, 1, 0]], idx = [];
    for (let j = 0; j <= seg; j++) { const a = 2 * Math.PI * j / seg; v.push([r * Math.cos(a), y + h, r * Math.sin(a)]); n.push([0, 1, 0]); }
    for (let j = 0; j < seg; j++) idx.push(0, j + 2, j + 1);
    const base = v.length;
    for (let j = 0; j <= seg; j++) { const a = 2 * Math.PI * j / seg, d = [Math.cos(a), 0, Math.sin(a)];
      v.push([r * d[0], y + h, r * d[2]]); n.push(d); v.push([r * d[0], y, r * d[2]]); n.push(d); }
    for (let j = 0; j < seg; j++) { const a = base + 2 * j; idx.push(a, a + 2, a + 1, a + 1, a + 2, a + 3); }
    return [v, n, idx];
  }
  function water(O, dir) {
    const d = norm(dir); let u = cross(d, [0, 1, 0]);
    if (Math.hypot(...u) < 1e-3) u = cross(d, [1, 0, 0]);
    u = norm(u); const a = 52 * Math.PI / 180;
    const Hs = [1, -1].map(s => add(O, sc(add(sc(d, Math.cos(a)), sc(u, s * Math.sin(a))), 0.05)));
    return [[['O', O], ...Hs.map(h => ['H', h])], Hs.map(h => [O, h])];
  }
  function glb(atoms, bonds, discColor) {
    const groups = {};
    const put = (k, g) => { const G = groups[k] || (groups[k] = [[], [], []]); const off = G[0].length;
      G[0].push(...g[0]); G[1].push(...g[1]); G[2].push(...g[2].map(i => i + off)); };
    atoms.forEach(([el, p]) => put(el, sphere(p, RAD[el])));
    bonds.forEach(([a, b]) => put('bond', cylinder(a, b, 0.009)));
    put('disc', disc(-0.27, 0.3));
    const colors = Object.assign({}, COL, {disc: discColor});
    const chunks = []; let len = 0; const views = [], accs = [], prims = [], mats = [];
    for (const [k, [v, n, idx]] of Object.entries(groups)) {
      const V = new Float32Array(v.flat()), N = new Float32Array(n.flat()), I = new Uint32Array(idx);
      const mn = [0, 1, 2].map(c => Math.min(...v.map(p => p[c]))), mx = [0, 1, 2].map(c => Math.max(...v.map(p => p[c])));
      for (const [arr, tgt] of [[V, 34962], [N, 34962], [I, 34963]]) {
        views.push({buffer: 0, byteOffset: len, byteLength: arr.byteLength, target: tgt});
        chunks.push(new Uint8Array(arr.buffer)); len += arr.byteLength; }
      const b = views.length - 3;
      accs.push({bufferView: b, componentType: 5126, count: v.length, type: 'VEC3', min: mn, max: mx},
                {bufferView: b + 1, componentType: 5126, count: n.length, type: 'VEC3'},
                {bufferView: b + 2, componentType: 5125, count: idx.length, type: 'SCALAR'});
      mats.push({name: k, pbrMetallicRoughness: {baseColorFactor: [...lin(colors[k]), 1], metallicFactor: 0.05, roughnessFactor: k === 'disc' ? 0.75 : 0.35}});
      prims.push({attributes: {POSITION: accs.length - 3, NORMAL: accs.length - 2}, indices: accs.length - 1, material: mats.length - 1});
    }
    const gltf = {asset: {version: '2.0', generator: 'hydrates-lab'}, scene: 0, scenes: [{nodes: [0]}], nodes: [{mesh: 0}],
      meshes: [{primitives: prims}], materials: mats, accessors: accs, bufferViews: views, buffers: [{byteLength: len}]};
    let js = new TextEncoder().encode(JSON.stringify(gltf)); const pad = (4 - js.length % 4) % 4;
    const jsP = new Uint8Array(js.length + pad).fill(0x20); jsP.set(js);
    const total = 12 + 8 + jsP.length + 8 + len, out = new Uint8Array(total), dv = new DataView(out.buffer);
    dv.setUint32(0, 0x46546C67, true); dv.setUint32(4, 2, true); dv.setUint32(8, total, true);
    dv.setUint32(12, jsP.length, true); dv.setUint32(16, 0x4E4F534A, true); out.set(jsP, 20);
    let o = 20 + jsP.length; dv.setUint32(o, len, true); dv.setUint32(o + 4, 0x004E4942, true); o += 8;
    for (const c of chunks) { out.set(c, o); o += c.length; }
    return URL.createObjectURL(new Blob([out], {type: 'model/gltf-binary'}));
  }
  function salt(metal, anions, anionBonds, waters) {
    const atoms = [[metal, [0, 0, 0]], ...anions], bonds = [...anionBonds];
    for (const [O, d, bonded] of waters) { const [a, b] = water(O, d); atoms.push(...a); bonds.push(...b); if (bonded) bonds.push([[0, 0, 0], O]); }
    return [atoms, bonds];
  }
  const R = 0.12;
  const cl = [['Cl', [0, R + 0.02, 0]], ['Cl', [0, -R - 0.02, 0]]], clb = cl.map(c => [[0, 0, 0], c[1]]);
  const w4 = [[[R, 0, 0], [1, 0, 0], 1], [[-R, 0, 0], [-1, 0, 0], 1], [[0, 0, R], [0, 0, 1], 1], [[0, 0, -R], [0, 0, -1], 1]];
  const w2 = [[[0.2, 0.12, 0.16], [1, 0.3, 1], 0], [[-0.2, -0.04, -0.16], [-1, 0.3, -1], 0]];
  const tet = [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]].map(t => sc(t, 1 / Math.sqrt(3)));
  const so4 = [['S', [0, 0.19, 0]], ...tet.map(t => ['O', add([0, 0.19, 0], sc(t, 0.055))])];
  const so4b = [...so4.slice(1).map(o => [[0, 0.19, 0], o[1]]), [[0, 0, 0], [0, 0.19, 0]]];
  const cuw = [...w4, [[0.18, -0.1, -0.16], [1, -0.3, -1], 0]];
  const bcl = [['Cl', [R + 0.03, 0, 0]], ['Cl', [-R - 0.03, 0, 0]]], bclb = bcl.map(c => [[0, 0, 0], c[1]]);
  const bw = [[[0, 0, R + 0.02], [0, 0, 1], 1], [[0, 0, -R - 0.02], [0, 0, -1], 1]];
  const defs = {
    CoCl2_6H2O: [salt('Co', cl, clb, [...w4, ...w2]), 'E06A9A'], CoCl2: [salt('Co', cl, clb, []), '3949AB'],
    CuSO4_5H2O: [salt('Cu', so4, so4b, cuw), '1E88E5'], CuSO4: [salt('Cu', so4, so4b, []), 'E8EAED'],
    BaCl2_2H2O: [salt('Ba', bcl, bclb, bw), 'D7E3EA'], BaCl2: [salt('Ba', bcl, bclb, []), 'F3F4F6']
  };
  const cache = {};
  window.saltModel = name => cache[name] || (cache[name] = glb(defs[name][0][0], defs[name][0][1], defs[name][1]));
})();
