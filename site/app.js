(() => {
'use strict';
const $ = s => document.querySelector(s);
const $$ = s => [...document.querySelectorAll(s)];
const SVGNS = 'http://www.w3.org/2000/svg';
const R = 8.314462618;           // J K-1 mol-1
const WIRE = 5.857;              // kJ g-1  (1400 cal g-1, course rubric value)
const BA_DU = 26.41;            // kJ g-1  (course rubric value)
const reduceMotion = window.matchMedia && matchMedia('(prefers-reduced-motion: reduce)').matches;

/* ---------- compounds ---------- */
const CMP = {
  ba: { name: 'Benzoic acid', formula: 'C₇H₆O₂', phase: 's', M: 122.12, du: -26.41, C: 7, H: 6, O: 2, O2: 7.5, CO2: 7, H2O: 3 },
  n:  { name: 'Naphthalene',  formula: 'C₁₀H₈',  phase: 's', M: 128.17, du: -40.19,  C: 10, H: 8, O: 0, O2: 12, CO2: 10, H2O: 4 },
  c:  { name: 'CCDT',         formula: 'C₁₂H₁₈', phase: 'l', M: 162.28, du: -45.10,  C: 12, H: 18, O: 0, O2: 16.5, CO2: 12, H2O: 9 }
};

/* ---------- storage ---------- */
const KEY = 'bombcal-v1';
let store = {};
try { store = JSON.parse(localStorage.getItem(KEY) || '{}') || {}; } catch (e) { store = {}; }
function persist() { try { localStorage.setItem(KEY, JSON.stringify(store)); } catch (e) { /* storage unavailable */ } }
function bindSaves(root = document) {
  root.querySelectorAll('[data-save]').forEach(elm => {
    if (elm.dataset.bound) return;
    elm.dataset.bound = '1';
    const id = elm.id;
    if (id in store) { if (elm.type === 'checkbox') elm.checked = !!store[id]; else elm.value = store[id]; }
    const ev = (elm.type === 'checkbox' || elm.tagName === 'SELECT') ? 'change' : 'input';
    elm.addEventListener(ev, () => { store[id] = elm.type === 'checkbox' ? elm.checked : elm.value; persist(); });
  });
}
function setVal(id, v) { const e = document.getElementById(id); if (!e) return; e.value = v; store[id] = v; persist(); e.dispatchEvent(new Event('input', { bubbles: true })); }

/* ---------- helpers ---------- */
const fmt = (x, d = 3) => Number.isFinite(x) ? x.toFixed(d) : '—';
const num = id => { const e = document.getElementById(id); if (!e || e.value.trim() === '') return NaN; return parseFloat(e.value); };
function svgEl(tag, attrs = {}, parent) { const e = document.createElementNS(SVGNS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; }
function svgText(parent, x, y, str, cls, anchor = 'start') { const t = svgEl('text', { x, y, class: cls, 'text-anchor': anchor }, parent); t.textContent = str; return t; }
function mulberry32(a) { return function () { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
function esc(s) { return String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c])); }
function parseFrac(s) {
  s = String(s || '').trim(); if (!s) return NaN;
  let m = s.match(/^(\d+)\s+(\d+)\/(\d+)$/); if (m) return +m[1] + (+m[2]) / (+m[3]);
  m = s.match(/^(-?\d+(?:\.\d+)?)\/(\d+(?:\.\d+)?)$/); if (m) return (+m[1]) / (+m[2]);
  return Number(s);
}
function fb(el, kind, html) { el.className = 'fb ' + (kind || ''); el.innerHTML = html || ''; }

/* ---------- router ---------- */
const VIEWS = ['overview', 'prelab', 'manual', 'postlab', 'references'];
function route() {
  let h = location.hash.slice(1);
  if (!VIEWS.includes(h)) h = 'overview';
  $$('.view').forEach(v => { v.hidden = v.dataset.view !== h; });
  $$('.tabs a').forEach(a => { if (a.dataset.tab === h) a.setAttribute('aria-current', 'page'); else a.removeAttribute('aria-current'); });
  window.scrollTo(0, 0);
}
window.addEventListener('hashchange', route);
document.addEventListener('click', e => {
  const a = e.target.closest('[data-jump]');
  if (!a) return;
  e.preventDefault();
  const t = document.getElementById(a.dataset.jump);
  if (t) t.scrollIntoView({ behavior: reduceMotion ? 'auto' : 'smooth', block: 'start' });
});

/* =====================================================================
   SIMULATION MODEL (isoperibol calorimeter)
   ===================================================================== */
function makeRun(key, m_mg, w_mg, o = {}) {
  const c = CMP[key];
  const Ccal = o.C || 2250;            // J K-1
  const Tj = o.Tj ?? 24.40;            // jacket, °C
  let T = o.T0 ?? 23.200;              // bucket, °C
  const Q = (m_mg / 1000) * Math.abs(c.du) * 1000 + (w_mg / 1000) * WIRE * 1000; // J
  const tau = o.tau || 22, k = 1.0e-4, stir = 2e-6;
  const rnd = mulberry32(o.seed || Math.floor(Math.random() * 1e9));
  const gauss = () => { let u = 0, v = 0; while (u === 0) u = rnd(); while (v === 0) v = rnd(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); };
  let t = 0, tf = null, next = 0;
  const data = { t: [], T: [] };
  function sample() { while (next <= t + 1e-9) { data.t.push(next); data.T.push(+(T + gauss() * 0.0004).toFixed(4)); next += 5; } }
  sample();
  function advance(dt) {
    const end = t + dt;
    while (t < end - 1e-9) {
      const h = Math.min(1, end - t);
      let src = 0;
      if (tf !== null && t + h > tf) { const x0 = Math.max(0, t - tf) / tau, x1 = (t + h - tf) / tau; src = Q / Ccal * (Math.exp(-x0) - Math.exp(-x1)); }
      T += src + (k * (Tj - T) + stir) * h;
      t += h;
      sample();
    }
  }
  return { advance, fire() { if (tf === null) tf = t; }, get t() { return t; }, get T() { return T; }, get fired() { return tf !== null; }, data, Q };
}
function exampleRun(key) {
  const p = { ba: [180.2, 6.0, 11], n: [101.5, 5.8, 23], c: [100.4, 6.1, 37] }[key];
  const r = makeRun(key, p[0], p[1], { seed: p[2] });
  r.advance(360); r.fire(); r.advance(840);
  return r.data;
}

/* ---------- regression & ΔT analysis ---------- */
function linfit(xs, ys) {
  const n = xs.length; let sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (let i = 0; i < n; i++) { sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; }
  const d = n * sxx - sx * sx; const b = d ? (n * sxy - sx * sy) / d : 0; const a = (sy - b * sx) / n;
  return { a, b, at: x => a + b * x };
}
function autoBounds(d) {
  const T = d.T, n = T.length; let iMax = 0;
  for (let i = 0; i < n; i++) if (T[i] > T[iMax]) iMax = i;
  const base = T[0], rise = T[iMax] - base;
  let iF = 0; for (let i = 0; i < n; i++) { if (T[i] > base + 0.03 * rise) { iF = Math.max(1, i - 1); break; } }
  const iP = Math.min(n - 3, iMax + Math.round(0.12 * (n - iMax)));
  return { iF, iP };
}
function analyze(d, iF, iP, frac) {
  const t = d.t, T = d.T, n = t.length;
  iF = Math.max(2, Math.min(iF, n - 4)); iP = Math.max(iF + 1, Math.min(iP, n - 3));
  const pre = linfit(t.slice(0, iF + 1), T.slice(0, iF + 1));
  const post = linfit(t.slice(iP), T.slice(iP));
  const tf = t[iF], tp = t[iP];
  const Ta = pre.at(tf), Tc = post.at(tp);
  const target = Ta + frac * (Tc - Ta);
  let tb = tf;
  for (let i = iF; i < iP; i++) {
    if (T[i] <= target && T[i + 1] >= target) { tb = t[i] + (target - T[i]) / ((T[i + 1] - T[i]) || 1e-9) * (t[i + 1] - t[i]); break; }
    if (i === iP - 1) tb = tf + frac * (tp - tf);
  }
  return { pre, post, tf, tp, tb, iF, iP, dT: post.at(tb) - pre.at(tb), raw: Tc - Ta };
}

/* ---------- chart ---------- */
function niceStep(span, count) { const s0 = span / count; const mag = Math.pow(10, Math.floor(Math.log10(s0))); const e = s0 / mag; return (e >= 7.5 ? 10 : e >= 3.5 ? 5 : e >= 1.5 ? 2 : 1) * mag; }
function ticks(min, max, count) { const st = niceStep(max - min, count); const out = []; for (let v = Math.ceil(min / st) * st; v <= max + st * 1e-6; v += st) out.push(+v.toFixed(10)); return { out, st }; }
const decs = st => Math.max(0, Math.ceil(-Math.log10(st) - 1e-9));
function plot(host, d, o = {}) {
  const W = 720, H = o.h || 340, ml = 70, mr = 20, mt = 18, mb = 48;
  host.innerHTML = '';
  const box = document.createElement('div'); box.className = 'scroll-x'; host.appendChild(box);
  const svg = svgEl('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': o.label || 'Temperature versus time', style: 'min-width:540px' }, box);
  if (!d || !d.t.length) { svgText(svg, W / 2, H / 2, o.empty || 'No data yet', 'ax-lbl', 'middle'); return; }
  const tm = d.t.map(x => x / 60);
  const xmin = o.xmin ?? tm[0], xmax = Math.max(o.xmax ?? tm[tm.length - 1], xmin + 0.5);
  let ymin = Math.min(...d.T), ymax = Math.max(...d.T);
  if (o.fit) { ymin = Math.min(ymin, o.fit.pre.at(o.fit.tb)); ymax = Math.max(ymax, o.fit.post.at(o.fit.tb)); }
  const pad = (ymax - ymin) * 0.08 || 0.01; ymin -= pad; ymax += pad;
  const X = v => ml + (v - xmin) / (xmax - xmin) * (W - ml - mr);
  const Y = v => H - mb - (v - ymin) / (ymax - ymin) * (H - mt - mb);
  if (o.regions) o.regions.forEach(r => {
    svgEl('rect', { class: 'reg', x: X(r.a), y: mt, width: Math.max(0, X(r.b) - X(r.a)), height: H - mt - mb }, svg);
    svgText(svg, (X(r.a) + X(r.b)) / 2, mt + 16, r.label, 'reg-t', 'middle');
  });
  const xt = ticks(xmin, xmax, 8), yt = ticks(ymin, ymax, 6);
  xt.out.forEach(v => { svgEl('line', { class: 'grid', x1: X(v), x2: X(v), y1: mt, y2: H - mb }, svg); svgText(svg, X(v), H - mb + 17, v.toFixed(decs(xt.st)), 'tick', 'middle'); });
  yt.out.forEach(v => { svgEl('line', { class: 'grid', x1: ml, x2: W - mr, y1: Y(v), y2: Y(v) }, svg); svgText(svg, ml - 8, Y(v) + 4, v.toFixed(decs(yt.st)), 'tick', 'end'); });
  svgEl('line', { class: 'axis', x1: ml, x2: ml, y1: mt, y2: H - mb }, svg);
  svgEl('line', { class: 'axis', x1: ml, x2: W - mr, y1: H - mb, y2: H - mb }, svg);
  svgText(svg, (ml + W - mr) / 2, H - 8, 'Time (min)', 'ax-lbl', 'middle');
  const yl = svgText(svg, 16, (mt + H - mb) / 2, 'Temperature (°C)', 'ax-lbl', 'middle'); yl.setAttribute('transform', `rotate(-90 16 ${(mt + H - mb) / 2})`);
  if (o.fit) {
    const f = o.fit, s = v => v / 60;
    const seg = (cls, fn, a, b) => svgEl('line', { class: cls, x1: X(s(a)), y1: Y(fn.at(a)), x2: X(s(b)), y2: Y(fn.at(b)) }, svg);
    seg('fitreg', f.pre, d.t[0], f.tf); seg('fit', f.pre, f.tf, f.tb);
    seg('fitreg', f.post, f.tp, d.t[d.t.length - 1]); seg('fit', f.post, f.tb, f.tp);
    svgEl('line', { class: 'vline', x1: X(s(f.tf)), x2: X(s(f.tf)), y1: mt, y2: H - mb }, svg);
    svgEl('line', { class: 'vline', x1: X(s(f.tp)), x2: X(s(f.tp)), y1: mt, y2: H - mb }, svg);
  }
  const pts = tm.map((v, i) => `${X(v).toFixed(1)},${Y(d.T[i]).toFixed(1)}`).join(' ');
  svgEl('polyline', { class: 'trace', points: pts }, svg);
  if (o.fit) {
    const f = o.fit, xb = X(f.tb / 60), y1 = Y(f.pre.at(f.tb)), y2 = Y(f.post.at(f.tb));
    svgEl('line', { class: 'dt', x1: xb, x2: xb, y1, y2 }, svg);
    svgEl('line', { class: 'dt', x1: xb - 6, x2: xb + 6, y1, y2: y1 }, svg);
    svgEl('line', { class: 'dt', x1: xb - 6, x2: xb + 6, y1: y2, y2 }, svg);
    svgText(svg, xb + 10, (y1 + y2) / 2, `ΔT = ${f.dT.toFixed(4)} K`, 'dt-t');
  }
  return { svg, X, Y };
}

/* =====================================================================
   OVERVIEW · interactive calorimeter
   ===================================================================== */
const PARTS = {
  bomb: ['1 · Oxygen bomb', 'A thick-walled stainless-steel vessel with a screw cap. Because it is sealed and rigid, combustion happens at constant volume, so the heat released equals ΔU, not ΔH. It is charged with pure O₂ (about 30 atm) so combustion is fast and complete.'],
  sample: ['2 · Sample cup and pellet', 'Solids are pressed into a pellet so they burn evenly and don\'t scatter. The liquid CCDT is pipetted into a platinum boat that is weighed before and after filling. Total heat released must stay under 1200 cal.'],
  fuse: ['3 · Electrodes and fuse wire', 'A short pulse of current heats a fine fuse wire that touches the sample and ignites it. Part of the wire burns and releases its own heat (1400 cal g⁻¹), so you weigh the unburned pieces and correct for it.'],
  valve: ['4 · Oxygen valve', 'Used to fill the bomb with oxygen before the run and to release the pressure slowly afterwards. Never unscrew the cap while the bomb is still pressurized.'],
  water: ['5 · Bucket and water', 'The bomb sits fully submerged in an exactly measured volume of water. The water, bucket and bomb together are the energy equivalent C. Use the same water volume in every run, because changing it changes C.'],
  jacket: ['6 · Insulated jacket', 'An air gap and an insulated jacket at constant temperature (isoperibol operation) slow heat exchange with the room. Some heat still leaks, which is why you correct ΔT by extrapolating the fore and after drifts.'],
  stirrer: ['7 · Stirrer', 'Keeps the water at a uniform temperature so the probe reads the true bucket temperature. Stirring adds a small, steady heat input that shows up in the drift lines.'],
  thermistor: ['8 · Thermistor probe', 'A semiconductor sensor whose resistance falls as temperature rises (negative temperature coefficient). The microprocessor converts resistance to temperature and logs it for later analysis.'],
  ignition: ['9 · Ignition leads', 'Carry the firing current from the ignition unit to the bomb electrodes. Attach them before closing the cover and stand clear when you fire.']
};
function initCal() {
  const info = $('#calInfo');
  $$('#calSvg .part').forEach(g => {
    const show = () => {
      $$('#calSvg .part').forEach(p => p.classList.remove('active'));
      g.classList.add('active');
      const [h, p] = PARTS[g.dataset.part];
      info.innerHTML = `<h4>${h}</h4><p style="margin:.3rem 0 0">${p}</p>`;
    };
    g.addEventListener('click', show);
    g.addEventListener('keydown', e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); show(); } });
  });
  const lcd = $('#calLcd'), flash = $('#calFlash'), pellet = $('#calPellet'), prop = $('#calProp'), line = $('#calSparkLine');
  const T0 = 23.2, dT = 2.13;
  let raf = null;
  const Tof = s => { // s: simulated seconds 0..900, fire at 120
    let T = T0 + 0.00012 * s;
    if (s > 120) { const x = s - 120; T += dT * (1 - Math.exp(-x / 30)) - 0.00009 * Math.max(0, x - 120); }
    return T;
  };
  function reset() {
    cancelAnimationFrame(raf); raf = null;
    lcd.innerHTML = T0.toFixed(3) + '<small>°C</small>';
    flash.setAttribute('opacity', 0); pellet.setAttribute('opacity', 1); line.setAttribute('points', '');
    prop.removeAttribute('transform'); $('#calFire').disabled = false;
  }
  function fire() {
    reset(); $('#calFire').disabled = true;
    const start = performance.now(), dur = reduceMotion ? 1 : 7000, pts = [];
    function frame(now) {
      const p = Math.min(1, (now - start) / dur), s = p * 900, T = Tof(s);
      lcd.innerHTML = T.toFixed(3) + '<small>°C</small>';
      pts.push(`${(4 + p * 152).toFixed(1)},${(40 - (T - T0) / 2.3 * 34).toFixed(1)}`);
      line.setAttribute('points', pts.join(' '));
      const fx = (s - 120) / 40;
      if (fx > 0 && fx < 1) { flash.setAttribute('opacity', (1 - fx).toFixed(2)); flash.setAttribute('r', (10 + fx * 36).toFixed(1)); }
      else flash.setAttribute('opacity', 0);
      if (s > 125) pellet.setAttribute('opacity', 0.15);
      if (!reduceMotion) prop.setAttribute('transform', `translate(410 0) scale(${Math.cos(now / 70).toFixed(2)} 1) translate(-410 0)`);
      if (p < 1) raf = requestAnimationFrame(frame); else { $('#calFire').disabled = false; prop.removeAttribute('transform'); }
    }
    raf = requestAnimationFrame(frame);
  }
  $('#calFire').addEventListener('click', fire);
  $('#calReset').addEventListener('click', reset);
}

/* =====================================================================
   PRE-LAB
   ===================================================================== */
function ring(cx, cy, r, n, rot) { const p = []; for (let i = 0; i < n; i++) { const a = (rot + i * 360 / n) * Math.PI / 180; p.push([cx + r * Math.cos(a), cy + r * Math.sin(a)]); } return p; }
function bond(svg, p, q, dbl, center) {
  svgEl('line', { class: 'mol-line', x1: p[0], y1: p[1], x2: q[0], y2: q[1] }, svg);
  if (dbl) {
    const mx = (p[0] + q[0]) / 2, my = (p[1] + q[1]) / 2; let dx = center[0] - mx, dy = center[1] - my; const L = Math.hypot(dx, dy); dx = dx / L * 5; dy = dy / L * 5;
    const sh = 0.18; const a = [p[0] + (q[0] - p[0]) * sh + dx, p[1] + (q[1] - p[1]) * sh + dy], b = [q[0] + (p[0] - q[0]) * sh + dx, q[1] + (p[1] - q[1]) * sh + dy];
    svgEl('line', { class: 'mol-line', x1: a[0], y1: a[1], x2: b[0], y2: b[1] }, svg);
  }
}
function drawRing(svg, pts, dblIdx, center) { pts.forEach((p, i) => bond(svg, p, pts[(i + 1) % pts.length], dblIdx.includes(i), center)); }
function molecules() {
  const host = $('#molCards');
  const defs = [
    ['Benzoic acid', 'C₇H₆O₂ · calibrant', s => {
      const c = [55, 55], p = ring(c[0], c[1], 26, 6, 0); drawRing(s, p, [1, 3, 5], c);
      const C = [p[0][0] + 24, p[0][1]]; svgEl('line', { class: 'mol-line', x1: p[0][0], y1: p[0][1], x2: C[0], y2: C[1] }, s);
      svgEl('line', { class: 'mol-line', x1: C[0], y1: C[1], x2: C[0] + 13, y2: C[1] - 19 }, s);
      svgEl('line', { class: 'mol-line', x1: C[0] + 5, y1: C[1] + 1, x2: C[0] + 17, y2: C[1] - 16 }, s);
      svgEl('line', { class: 'mol-line', x1: C[0], y1: C[1], x2: C[0] + 13, y2: C[1] + 19 }, s);
      svgText(s, C[0] + 19, C[1] - 26, 'O', 'mol-txt', 'middle'); svgText(s, C[0] + 24, C[1] + 26, 'OH', 'mol-txt', 'middle');
    }],
    ['Naphthalene', 'C₁₀H₈ · solid unknown', s => {
      const r = 25, L = [75 - r * Math.sqrt(3) / 2, 55], Rr = [75 + r * Math.sqrt(3) / 2, 55];
      const A = a => [L[0] + r * Math.cos(a * Math.PI / 180), L[1] + r * Math.sin(a * Math.PI / 180)], B = a => [Rr[0] + r * Math.cos(a * Math.PI / 180), Rr[1] + r * Math.sin(a * Math.PI / 180)];
      const s1 = A(-30), s2 = A(30), a1 = A(-90), a2 = A(210), a3 = A(150), a4 = A(90), b1 = B(-90), b2 = B(-30), b3 = B(30), b4 = B(90);
      bond(s, s1, s2, true, L); bond(s, s1, a1, false, L); bond(s, a1, a2, true, L); bond(s, a2, a3, false, L); bond(s, a3, a4, true, L); bond(s, a4, s2, false, L);
      bond(s, s1, b1, false, Rr); bond(s, b1, b2, true, Rr); bond(s, b2, b3, false, Rr); bond(s, b3, b4, true, Rr); bond(s, b4, s2, false, Rr);
    }],
    ['CCDT', 'C₁₂H₁₈ · liquid unknown', s => { const c = [75, 55], p = ring(c[0], c[1], 42, 12, -90); drawRing(s, p, [0, 4, 8], c); }],
    ['Benzene', 'C₆H₆ · literature', s => { const c = [75, 55], p = ring(c[0], c[1], 32, 6, -90); drawRing(s, p, [0, 2, 4], c); }],
    ['Cyclohexane', 'C₆H₁₂ · literature', s => { const c = [75, 55], p = ring(c[0], c[1], 32, 6, -90); drawRing(s, p, [], c); }]
  ];
  defs.forEach(([name, sub, draw]) => {
    const d = document.createElement('div'); d.className = 'mol';
    const s = svgEl('svg', { viewBox: '0 0 150 110', role: 'img', 'aria-label': name + ' skeletal structure' }); draw(s);
    d.appendChild(s); d.insertAdjacentHTML('beforeend', `<b>${name}</b><span>${sub}</span>`); host.appendChild(d);
  });
  host.insertAdjacentHTML('afterend', '<p class="saved-note">CCDT is drawn as a flat 12-membered ring with double bonds at the 1, 5 and 9 positions. The real molecule is puckered, with two <i>trans</i> and one <i>cis</i> double bond.</p>');
}
function balancer() {
  const host = $('#balancer');
  Object.entries(CMP).forEach(([k, c]) => {
    const row = document.createElement('div'); row.className = 'rxn';
    row.innerHTML = `<span class="name">${c.name}</span>
      <span>1 ${c.formula}(${c.phase}) +</span><input type="text" inputmode="decimal" id="bal_${k}_O2" aria-label="${c.name} O2 coefficient" data-save><span>O₂(g) →</span>
      <input type="text" inputmode="decimal" id="bal_${k}_CO2" aria-label="${c.name} CO2 coefficient" data-save><span>CO₂(g) +</span>
      <input type="text" inputmode="decimal" id="bal_${k}_H2O" aria-label="${c.name} H2O coefficient" data-save><span>H₂O</span>
      <select id="bal_${k}_ph" aria-label="${c.name} water phase" data-save><option value="">(phase)</option><option value="l">(l)</option><option value="g">(g)</option></select>
      <span class="mark" id="bal_${k}_atoms"></span>`;
    host.appendChild(row);
  });
  host.insertAdjacentHTML('beforeend', '<div class="btn-row"><button id="balCheck">Check all</button><span class="saved-note">The atom tally updates as you type.</span></div><div class="fb" id="balFb"></div>');
  const tally = () => Object.entries(CMP).forEach(([k, c]) => {
    const o2 = parseFrac($(`#bal_${k}_O2`).value), co2 = parseFrac($(`#bal_${k}_CO2`).value), h2o = parseFrac($(`#bal_${k}_H2O`).value);
    const out = $(`#bal_${k}_atoms`);
    if ([o2, co2, h2o].some(x => !Number.isFinite(x))) { out.textContent = ''; return; }
    const L = { C: c.C, H: c.H, O: c.O + 2 * o2 }, Rt = { C: co2, H: 2 * h2o, O: 2 * co2 + h2o };
    out.innerHTML = ['C', 'H', 'O'].map(a => { const ok = Math.abs(L[a] - Rt[a]) < 1e-9; return `<span style="color:var(--${ok ? 'ok' : 'bad'})">${a} ${+L[a].toFixed(2)}|${+Rt[a].toFixed(2)}</span>`; }).join(' · ');
  });
  host.addEventListener('input', tally); host.addEventListener('change', tally);
  setTimeout(tally, 0);
  $('#balCheck').addEventListener('click', () => {
    const msgs = []; let all = true;
    Object.entries(CMP).forEach(([k, c]) => {
      const o2 = parseFrac($(`#bal_${k}_O2`).value), co2 = parseFrac($(`#bal_${k}_CO2`).value), h2o = parseFrac($(`#bal_${k}_H2O`).value), ph = $(`#bal_${k}_ph`).value;
      if (Math.abs(o2 - c.O2) < 1e-9 && Math.abs(co2 - c.CO2) < 1e-9 && Math.abs(h2o - c.H2O) < 1e-9) {
        if (ph === 'l') msgs.push(`✓ ${c.name}: balanced, with the correct water phase.`);
        else { all = false; msgs.push(`◐ ${c.name}: balanced. Now set the water phase. At 25 °C and 1 bar, water's standard state is liquid.`); }
      } else { all = false; msgs.push(`✗ ${c.name}: not balanced yet. Balance C first, then H, then O, and remember that each O₂ supplies two O atoms.`); }
    });
    fb($('#balFb'), all ? 'ok' : 'warn', msgs.join('<br>'));
  });
}
function dnWidget() {
  const sel = $('#dnCmp');
  Object.entries(CMP).forEach(([k, c]) => sel.insertAdjacentHTML('beforeend', `<option value="${k}">${c.name} (${c.formula})</option>`));
  const upd = () => {
    const c = CMP[sel.value], T = parseFloat($('#dnT').value); $('#dnTv').textContent = T.toFixed(2);
    const dn = c.CO2 - c.O2, U = c.du * c.M, corr = dn * R * T / 1000, H = U + corr;
    $('#dnOut').innerHTML = `Δn<sub>gas</sub> = ${c.CO2} − ${c.O2} = <b>${dn}</b> mol &nbsp;·&nbsp; Δn<sub>gas</sub>RT = <b>${corr.toFixed(3)}</b> kJ mol⁻¹<br>ΔcU ≈ ${U.toFixed(1)} kJ mol⁻¹ &nbsp;→&nbsp; ΔcH ≈ <b>${H.toFixed(1)}</b> kJ mol⁻¹`;
    $('#dnNote').textContent = `The correction is ${(Math.abs(corr / U) * 100).toFixed(3)} % of ΔcU. That is small, but for CCDT (${Math.abs(CMP.c.CO2 - CMP.c.O2)} mol of gas lost per mole burned) it is about as large as the uncertainties you will propagate in Part 5. Water counts as liquid, so it doesn't enter Δn.`;
  };
  sel.addEventListener('change', upd); $('#dnT').addEventListener('input', upd); upd();
}
function calibWidget() {
  const Ctrue = 2215;
  const upd = () => {
    const m = parseFloat($('#cbM').value), w = parseFloat($('#cbW').value);
    $('#cbMv').textContent = m.toFixed(4); $('#cbWv').textContent = w.toFixed(1);
    const qBA = m * BA_DU * 1000, qW = w / 1000 * WIRE * 1000, dT = (qBA + qW) / Ctrue;
    $('#cbOut').innerHTML = `q from benzoic acid = ${m.toFixed(4)} g × 26.41 kJ g⁻¹ = <b>${qBA.toFixed(1)}</b> J<br>q from fuse wire = ${(w / 1000).toFixed(4)} g × 5.857 kJ g⁻¹ = <b>${qW.toFixed(1)}</b> J<br>Measured adiabatic ΔT = <b>${dT.toFixed(4)}</b> K`;
    calibWidget.state = { qBA, qW, dT: +dT.toFixed(4) };
  };
  ['cbM', 'cbW'].forEach(id => $('#' + id).addEventListener('input', () => { upd(); fb($('#cbFb'), '', ''); }));
  upd();
  $('#cbCheck').addEventListener('click', () => {
    const a = parseFloat($('#cbAns').value), s = calibWidget.state, ref = (s.qBA + s.qW) / s.dT;
    if (!Number.isFinite(a)) return fb($('#cbFb'), 'warn', 'Enter your value of C first.');
    if (Math.abs(a - ref) / ref < 0.004) return fb($('#cbFb'), 'ok', `✓ Correct, C ≈ ${ref.toFixed(0)} J K⁻¹. Change the pellet mass: ΔT changes, but C stays the same because it belongs to the calorimeter, not the sample.`);
    if (Math.abs(a - s.qBA / s.dT) / ref < 0.004) return fb($('#cbFb'), 'warn', 'Close. You left out the heat from the fuse wire, which makes C slightly too small.');
    if (Math.abs(a * 1000 - ref) / ref < 0.01) return fb($('#cbFb'), 'warn', 'Check your units. It looks like you reported kJ K⁻¹ instead of J K⁻¹.');
    fb($('#cbFb'), 'bad', 'Not quite. Add both heat terms and divide by ΔT.');
  });
}
function hess() {
  const T = [
    'Press <b>Next step</b> to walk through the cycle one arrow at a time.',
    '<b>Step 1. Get CCDT into the gas phase.</b> You burn <i>liquid</i> CCDT, but the bond-counting argument holds for gas-phase molecules. Raising the liquid to the gas costs ΔvapH(CCDT), which you estimate with Trouton\'s rule.',
    '<b>Step 2. Your measurement.</b> The bomb gives ΔcU for liquid CCDT. Convert it to ΔcH with Δn<sub>gas</sub>RT. Condensing the gas (−ΔvapH) and then burning the liquid (ΔcH(l)) is the same as burning the gas directly: ΔcH(CCDT, g) = ΔcH(l) − ΔvapH.',
    '<b>Step 3. The literature path.</b> Because C₆H₆ + C₆H₁₂ = C₁₂H₁₈, burning one mole each of gaseous benzene and cyclohexane gives the same products. Use the liquid values from NIST and correct each with its own ΔvapH.',
    '<b>Step 4. The gap is the resonance energy.</b> Both paths end at 12 CO₂ + 9 H₂O, so the difference between the starting levels is E<sub>res</sub>. Real benzene sits <i>lower</i> than "cyclohexatriene", so its combustion releases less heat and E<sub>res</sub> comes out positive.'
  ];
  let k = 0;
  const show = () => {
    $$('#hess .hstep').forEach(g => { const s = +g.dataset.s; g.classList.toggle('on', s === k); g.classList.toggle('dim', k > 0 && s !== k); });
    $('#hText').innerHTML = `<p>${T[k]}</p>`; $('#hCount').textContent = `Step ${k} of 4`;
    $('#hPrev').disabled = k === 0; $('#hNext').disabled = k === 4;
  };
  $('#hNext').addEventListener('click', () => { if (k < 4) { k++; show(); } });
  $('#hPrev').addEventListener('click', () => { if (k > 0) { k--; show(); } });
  show();
}
function sorter() {
  const items = [
    ['You forget the fuse-wire correction on the naphthalene run.', 'sys', 'over', 'The heat from the burning wire gets counted as heat from the sample.'],
    ['Soot is left in the cup after the CCDT run (incomplete combustion).', 'sys', 'under', 'Less heat is released than complete combustion of the weighed mass would give, so |ΔcU| comes out too small.'],
    ['The thermistor reading scatters by ±0.0005 °C from point to point.', 'rand', 'either', 'Noise scatters both ways and averages down when you fit many points.'],
    ['The bucket holds 2 mL more water in the naphthalene run than in the calibration run.', 'sys', 'under', 'The true energy equivalent is larger than your calibrated C, so the same heat gives a smaller ΔT.'],
    ['Air left in the bomb lets some N₂ form nitric acid, and you don\'t correct for it.', 'sys', 'over', 'Forming HNO₃ is exothermic, and that extra heat is credited to the sample. Parr describes an acid titration correction for this.'],
    ['The balance reads 0.3 mg high when you weigh the naphthalene pellet.', 'sys', 'under', 'Dividing by a mass that is too large makes the energy per gram too small, by about 0.3 % here.'],
    ['Lab partners draw the extrapolation lines slightly differently.', 'rand', 'either', 'This is judgment scatter with no preferred direction. A consistent fitting rule reduces it.']
  ];
  const host = $('#sorter');
  items.forEach((it, i) => {
    const d = document.createElement('div'); d.className = 'sort-item';
    d.innerHTML = `<p>${it[0]}</p>
      <div class="choice-row" data-g="type" role="group" aria-label="Error type"><button data-v="sys" aria-pressed="false">Systematic</button><button data-v="rand" aria-pressed="false">Random</button></div>
      <div class="choice-row" data-g="dir" role="group" aria-label="Effect on |ΔcU|"><button data-v="over" aria-pressed="false">|ΔcU| too large</button><button data-v="under" aria-pressed="false">|ΔcU| too small</button><button data-v="either" aria-pressed="false">Either way</button></div>
      <div class="fb"></div>`;
    const pick = {};
    d.addEventListener('click', e => {
      const b = e.target.closest('button'); if (!b) return;
      const g = b.parentElement.dataset.g; pick[g] = b.dataset.v;
      b.parentElement.querySelectorAll('button').forEach(x => x.setAttribute('aria-pressed', x === b ? 'true' : 'false'));
      if (pick.type && pick.dir) {
        const ok = pick.type === it[1] && pick.dir === it[2];
        fb(d.querySelector('.fb'), ok ? 'ok' : 'bad', (ok ? '✓ ' : '✗ Not quite. ') + `${it[1] === 'sys' ? 'Systematic' : 'Random'}, ${{ over: '|ΔcU| too large', under: '|ΔcU| too small', either: 'either direction' }[it[2]]}. ${it[3]}`);
      }
    });
    host.appendChild(d);
  });
}
function quiz() {
  const Q = [
    ['In a sealed, rigid bomb, the heat released by the reaction equals…', ['ΔH of combustion', 'ΔU of combustion', 'ΔG of combustion', 'TΔS of combustion'], 1, 'At constant volume no PV work is done, so q<sub>V</sub> = ΔU.'],
    ['Why is benzoic acid used for calibration?', ['It burns at a lower temperature', 'Its specific energy of combustion is certified, and it is a stable, non-hygroscopic solid that is easy to pellet', 'It doesn\'t need oxygen', 'Its Δn<sub>gas</sub> is zero'], 1, 'NIST certifies its energy of combustion (SRM 39j), and it weighs and burns reproducibly.'],
    ['For C₁₀H₈(s) + 12 O₂(g) → 10 CO₂(g) + 4 H₂O(l), Δn<sub>gas</sub> is…', ['−2', '+2', '−6', '+4'], 0, '10 mol of gaseous CO₂ minus 12 mol of O₂ gives −2. Liquid water is not counted.'],
    ['ΔcU(naphthalene) = −5151 kJ mol⁻¹ at 298 K. ΔcH is about…', ['−5146 kJ mol⁻¹', '−5151 kJ mol⁻¹', '−5156 kJ mol⁻¹', '−5400 kJ mol⁻¹'], 2, 'Δn<sub>gas</sub>RT = (−2)(8.314)(298) J ≈ −4.96 kJ, so ΔH is about 5 kJ more negative than ΔU.'],
    ['As the temperature of a thermistor rises, its resistance…', ['increases', 'decreases', 'stays the same', 'drops to zero at a set point'], 1, 'Thermistors have a negative temperature coefficient of resistance.'],
    ['After ignition, why does the temperature slowly <i>fall</i> in the after period?', ['The sample keeps absorbing heat', 'The bucket is now warmer than the jacket and loses heat to it', 'The stirrer cools the water', 'The oxygen expands'], 1, 'Heat leaks from the warmer bucket to the cooler jacket. The extrapolation corrects for this.'],
    ['Trouton\'s rule estimates ΔvapH from…', ['the critical temperature', 'the normal boiling point, using ΔvapS ≈ 85–88 J K⁻¹ mol⁻¹', 'the melting point', 'the vapor density'], 1, 'ΔvapH ≈ T<sub>b</sub>·ΔvapS. For CCDT, (85 J K⁻¹ mol⁻¹)(504 K) ≈ 43 kJ mol⁻¹.']
  ];
  const host = $('#quiz'); let score = 0, answered = 0;
  host.innerHTML = '<div class="panel-title"><h4>Self-check</h4><span class="tag" id="quizScore">0 / ' + Q.length + '</span></div>';
  Q.forEach(([q, opts, ans, why]) => {
    const d = document.createElement('div'); d.className = 'mcq';
    d.innerHTML = `<p>${q}</p><div class="opts">${opts.map((o, i) => `<button data-i="${i}">${o}</button>`).join('')}</div><div class="fb"></div>`;
    let done = false;
    d.querySelector('.opts').addEventListener('click', e => {
      const b = e.target.closest('button'); if (!b || done) return; done = true; answered++;
      const i = +b.dataset.i; const btns = d.querySelectorAll('.opts button');
      btns[ans].classList.add('right'); if (i !== ans) b.classList.add('wrong'); else score++;
      fb(d.querySelector('.fb'), i === ans ? 'ok' : 'bad', (i === ans ? '✓ ' : '✗ ') + why);
      $('#quizScore').textContent = `${score} / ${Q.length}`;
    });
    host.appendChild(d);
  });
}
const PRELAB_QS = [
  ['1', 4, 'Write chemical equations for the combustion of (a) benzoic acid (1 pt), (b) naphthalene (1 pt), and (c) 1,5,9-<i>trans,trans,cis</i>-cyclododecatriene (CCDT) (1 pt). Balance each and include the phase of every component (1 pt).', 'Practice: Section B'],
  ['2', 4, '(a, 2 pts) Define the enthalpy of combustion (ΔcH) and the enthalpy of formation (ΔfH). (b, 2 pts) Explain how to find ΔfH for a compound from its ΔcH.', 'Think about: one mole, standard conditions and standard states, and Hess\'s law with tabulated values'],
  ['3', 3, 'Explain why we run a measurement with benzoic acid. What information does it give, and how is that information obtained quantitatively? Use equations in your answer.', 'Practice: Section D'],
  ['4', 4, '(a, 2 pts) Why is combustion carried out at constant volume rather than constant pressure in this experiment? (b, 2 pts) Explain how to convert energies of combustion from bomb calorimetry into enthalpies of combustion.', 'Practice: Sections A and C'],
  ['5', 3, 'Using Hess\'s law, explain how to find the resonance stabilization energy of benzene from ΔcH for benzene (C₆H₆), cyclohexane (C₆H₁₂) and CCDT (C₁₂H₁₈).', 'Practice: Section E (the cycle)'],
  ['6', 2, 'Describe a relevant <b>systematic</b> source of error for this experiment that is not a gross or human error (1 pt), and explain specifically how it would affect your results (1 pt).', 'Practice: the error sorter in Section E. Name the direction <i>and</i> follow it through to ΔcH.']
];
function prelabQs() {
  const host = $('#prelabQs');
  host.innerHTML = `<div class="panel"><div class="grid2"><div class="field"><label for="pl_name">Name</label><input type="text" id="pl_name" data-save autocomplete="name"></div><div class="field"><label for="pl_sec">Lab section / date</label><input type="text" id="pl_sec" data-save></div></div></div>` +
    PRELAB_QS.map(([n, p, q, h]) => `<div class="q"><div class="q-head"><h4>Question ${n}</h4><span class="pts">${p} pts</span></div><p>${q}</p><label for="pl_q${n}" class="sr">Answer to question ${n}</label><textarea id="pl_q${n}" data-save placeholder="Your answer…"></textarea><p class="saved-note">${h}</p></div>`).join('');
}

/* =====================================================================
   MANUAL
   ===================================================================== */
function chargeCheck() {
  const sel = $('#ccCmp');
  Object.entries(CMP).forEach(([k, c]) => sel.insertAdjacentHTML('beforeend', `<option value="${k}">${c.name} (≈ ${Math.abs(c.du).toFixed(1)} kJ g⁻¹)</option>`));
  const defaults = { ba: 180, n: 100, c: 100 };
  sel.addEventListener('change', () => { $('#ccM').value = defaults[sel.value]; upd(); });
  const upd = () => {
    const c = CMP[sel.value], m = num('ccM') / 1000, C = num('ccC');
    if (!(m > 0) || !(C > 0)) { $('#ccOut').textContent = 'Enter a mass and a heat capacity.'; return fb($('#ccFb'), '', ''); }
    const q = m * Math.abs(c.du) * 1000 + 0.006 * WIRE * 1000, cal = q / 4.184;
    $('#ccOut').innerHTML = `Heat released ≈ <b>${cal.toFixed(0)}</b> cal (${(q / 1000).toFixed(2)} kJ, including ~6 mg of fuse wire) &nbsp;·&nbsp; predicted ΔT ≈ <b>${(q / C).toFixed(2)}</b> K`;
    const maxMg = (1200 * 4.184 - 0.006 * WIRE * 1000) / (Math.abs(c.du) * 1000) * 1000;
    if (cal > 1200) fb($('#ccFb'), 'bad', `Over the 1200 cal limit. Use no more than about ${Math.floor(maxMg)} mg of ${c.name}.`);
    else if (cal > 1140) fb($('#ccFb'), 'warn', `Within 5 % of the limit. Fine if your mass is accurate, but ${Math.floor(maxMg * 0.95)} mg or less gives a safety margin.`);
    else fb($('#ccFb'), 'ok', `Safe charge. The maximum for ${c.name} is about ${Math.floor(maxMg)} mg.`);
  };
  ['ccM', 'ccC'].forEach(id => $('#' + id).addEventListener('input', upd)); upd();
}
function thermoFig() {
  const d = exampleRun('ba'), b = autoBounds(d), f = analyze(d, b.iF, b.iP, 0.6);
  plot($('#thermoFig'), d, { fit: f, label: 'Example thermogram for benzoic acid with extrapolated fore and after periods', regions: [
    { a: d.t[0] / 60, b: f.tf / 60, label: '① fore period' }, { a: f.tf / 60, b: f.tp / 60, label: '② rise' }, { a: f.tp / 60, b: d.t[d.t.length - 1] / 60, label: '③ after period' }] });
  $('#thermoFig').insertAdjacentHTML('beforeend', `<figcaption>Figure 1. Simulated thermogram for 0.1802 g of benzoic acid in a calorimeter with C ≈ 2250 J K⁻¹. Solid orange lines are least-squares fits to the fore and after periods. Dashed lines extrapolate them to the 60 % point (t = ${(f.tb / 60).toFixed(2)} min), where the vertical bar gives the adiabatic ΔT = ${f.dT.toFixed(4)} K. The raw rise without drift correction would be ${f.raw.toFixed(4)} K.</figcaption>`);
}
let virtualRun = null;
try { if (store._virtual) virtualRun = store._virtual; } catch (e) { /* ignore */ }
function simulator() {
  const sel = $('#simCmp');
  Object.entries(CMP).forEach(([k, c]) => sel.insertAdjacentHTML('beforeend', `<option value="${k}">${c.name}</option>`));
  const defaults = { ba: 180, n: 100, c: 100 };
  sel.addEventListener('change', () => { $('#simM').value = defaults[sel.value]; });
  let run = null, timer = null, last = 0, lastDraw = 0;
  const chart = $('#simChart');
  plot(chart, null, { empty: 'Press "Start recording" to begin a run' });
  const draw = () => plot(chart, run.data, { xmin: 0, xmax: Math.max(15, run.t / 60 + 1), label: 'Live virtual thermogram' });
  const setBtns = s => { $('#simStart').disabled = s !== 'idle'; $('#simFire').disabled = s !== 'rec'; $('#simStop').disabled = s === 'idle'; $('#simSkip').disabled = s === 'idle'; [sel, $('#simM'), $('#simW')].forEach(e => e.disabled = s !== 'idle'); };
  function tick(now) {
    const dt = Math.min(0.2, (now - last) / 1000); last = now;
    run.advance(dt * 30);
    $('#simLcd').innerHTML = run.T.toFixed(3) + '<small>°C</small>';
    $('#simClock').textContent = `t = ${Math.floor(run.t / 60)} min ${String(Math.floor(run.t % 60)).padStart(2, '0')} s${run.fired ? ' · fired' : ''}`;
    if (now - lastDraw > 150) { draw(); lastDraw = now; }
    if (run.t >= 1800) return stop('Run stopped automatically at 30 min.');
    timer = requestAnimationFrame(tick);
  }
  $('#simStart').addEventListener('click', () => {
    const k = sel.value, m = num('simM'), w = num('simW') || 0;
    if (!(m > 0)) return fb($('#simFb'), 'warn', 'Enter a sample mass.');
    const cal = ((m / 1000) * Math.abs(CMP[k].du) * 1000 + (w / 1000) * WIRE * 1000) / 4.184;
    if (cal > 1200) return fb($('#simFb'), 'bad', `This charge would release about ${cal.toFixed(0)} cal, which is over the 1200 cal limit. Your TA would not let you fire it. Reduce the mass.`);
    run = makeRun(k, m, w, { T0: 23.15 + Math.random() * 0.1 }); run.key = k; run.m = m; run.w = w;
    fb($('#simFb'), '', 'Recording the fore period. Wait for a steady, linear drift (at least 5 min of lab time) before you fire.'); $('#simFb').className = 'fb';
    setBtns('rec'); last = performance.now(); timer = requestAnimationFrame(tick);
  });
  $('#simFire').addEventListener('click', () => {
    if (!run) return;
    if (run.t < 180) fb($('#simFb'), 'warn', 'You fired after less than 3 min. A short fore period makes the drift line, and so ΔT, less certain.');
    else fb($('#simFb'), 'ok', 'Fired. Keep recording until the after period is clearly linear, about 8–10 min after the peak.');
    run.fire(); setBtns('fired');
  });
  $('#simSkip').addEventListener('click', () => { if (run) { run.advance(120); draw(); } });
  function stop(msg) {
    cancelAnimationFrame(timer); timer = null; setBtns('idle');
    if (!run) return; draw();
    const peakIdx = run.data.T.indexOf(Math.max(...run.data.T));
    const enough = run.fired && run.data.t.length - peakIdx > 60;
    virtualRun = { t: run.data.t.slice(), T: run.data.T.slice(), key: run.key, m: run.m, w: run.w };
    store._virtual = virtualRun; persist();
    const opt = $('#azSrc option[value="virtual"]'); if (opt) opt.disabled = false;
    fb($('#simFb'), enough ? 'ok' : 'warn', (msg ? msg + ' ' : '') + (run.fired ? (enough ? 'Run saved. In the post-lab analyzer, choose "My virtual run" to find ΔT.' : 'Run saved, but the after period is short. The after-period drift will be poorly defined.') : 'Run saved, but you never fired. There is no rise to analyze.'));
  }
  $('#simStop').addEventListener('click', () => stop());
}
const PROC = [
  ['Before you start', [
    'Put on safety glasses. Read the Parr operating pages and watch the procedure video on the course site.',
    'Check that the jacket temperature is stable and record it.',
    'Check that the oxygen tank regulator and filling connection are clean and free of oil.'
  ]],
  ['Prepare the sample', [
    'Solids: press a pellet of benzoic acid (≈ 0.18 g, never above 0.19 g) or naphthalene (≈ 0.10 g). Tap off loose fragments.',
    'Weigh the pellet on the analytical balance to ±0.0001 g and record it.',
    'Liquid CCDT: weigh the empty platinum boat, pipette in ≈ 0.10 g (≤ 0.11 g), and reweigh right away.',
    'Run the charge check to confirm the sample will release less than 1200 cal.'
  ]],
  ['Load the bomb', [
    'Cut a measured length of fuse wire, weigh it, and record the mass.',
    'Attach the wire between the two electrodes so it touches the pellet or dips into the liquid but does not touch the cup walls.',
    'Place the cup in its holder, check that the sealing ring and threads are clean, and screw the cap on hand-tight.'
  ]],
  ['Charge with oxygen', [
    'Connect the filling hose and admit O₂ slowly to the pressure in the Parr instructions (typically about 30 atm). Never exceed it.',
    'Close the valve and disconnect. If you suspect a leak, stop and call your TA.'
  ]],
  ['Assemble the calorimeter', [
    'Fill the bucket with the specified volume of water, measured exactly and the same for every run.',
    'Lower the bomb into the bucket, check for bubbles (a leak means do not fire), and attach the ignition leads.',
    'Close the cover, lower the thermistor, and start the stirrer.'
  ]],
  ['Record and fire', [
    'Let the system equilibrate (≈ 5 min), then start data logging.',
    'Record a fore period of at least 5 min with a steady, linear drift.',
    'Note the time and fire the bomb. Stand clear of the calorimeter.',
    'Record through the rise and keep recording for at least 8–10 min of linear after period.'
  ]],
  ['Disassemble and inspect', [
    'Stop the stirrer, raise the cover, and lift out the bomb.',
    'Release the pressure slowly through the valve before opening the cap.',
    'Inspect the cup. Soot or unburned sample means incomplete combustion, so discard the run.',
    'Collect and weigh the unburned fuse wire. Wire burned = initial − unburned.',
    'Rinse and dry the bomb interior and cup.'
  ]],
  ['Repeat', [
    'Run 2: naphthalene. Run 3: CCDT. Use the same water volume each time.',
    'Save or export each data-logger file and write the file name on the data sheet.'
  ]]
];
function procedure() {
  const host = $('#procedure'); let n = 0;
  host.innerHTML = PROC.map(([title, steps], pi) => `<div class="phase"><h3><span class="eyebrow">Phase ${pi + 1}</span>${title}</h3><ul class="steps">${steps.map(s => { const id = 'proc_' + (n++); return `<li><input type="checkbox" id="${id}" data-save><label for="${id}">${s}</label></li>`; }).join('')}</ul></div>`).join('');
  const total = n;
  const upd = () => { const boxes = host.querySelectorAll('input[type=checkbox]'); let c = 0; boxes.forEach(b => { b.closest('li').classList.toggle('done', b.checked); if (b.checked) c++; }); $('#procBar').style.width = (c / total * 100) + '%'; $('#procCount').textContent = `${c} of ${total} steps done`; };
  host.addEventListener('change', upd);
  $('#procReset').addEventListener('click', () => { host.querySelectorAll('input[type=checkbox]').forEach(b => { b.checked = false; store[b.id] = false; }); persist(); upd(); });
  procedure.upd = upd;
}
const DS_ROWS = [
  ['mass', 'Sample mass (g)'], ['boat', 'Pt boat, empty (g), CCDT only'], ['wi', 'Fuse wire, initial (g)'], ['wu', 'Fuse wire, unburned (g)'], ['wb', 'Fuse wire burned (g)', true],
  ['vol', 'Bucket water volume (mL)'], ['tj', 'Jacket temperature (°C)'], ['t0', 'Initial bucket temperature (°C)'], ['tf', 'Firing time (min:s)'], ['file', 'Data file name'], ['obs', 'Observations (soot? residue?)']
];
const DS_COLS = [['ba', 'Run 1 · Benzoic acid'], ['n', 'Run 2 · Naphthalene'], ['c', 'Run 3 · CCDT']];
function dataSheet() {
  const t = $('#dataSheet');
  t.innerHTML = `<thead><tr><th>Quantity</th>${DS_COLS.map(c => `<th>${c[1]}</th>`).join('')}</tr></thead><tbody>${DS_ROWS.map(([k, lbl, ro]) => `<tr><td>${lbl}</td>${DS_COLS.map(([c]) => `<td><input type="text" id="ds_${k}_${c}" aria-label="${lbl}, ${c}" ${ro ? 'readonly tabindex="-1" style="background:var(--sunk)"' : 'data-save'} ${['mass', 'boat', 'wi', 'wu', 'vol', 'tj', 't0'].includes(k) ? 'inputmode="decimal"' : ''}></td>`).join('')}</tr>`).join('')}</tbody>`;
  t.closest('.scroll-x').style.cssText = 'overflow-x:auto';
  t.querySelectorAll('input').forEach(i => { i.style.minWidth = '130px'; });
  const upd = () => DS_COLS.forEach(([c]) => { const wi = num(`ds_wi_${c}`), wu = num(`ds_wu_${c}`); $(`#ds_wb_${c}`).value = Number.isFinite(wi) && Number.isFinite(wu) ? (wi - wu).toFixed(4) : ''; });
  t.addEventListener('input', upd); dataSheet.upd = upd;
  const row = $('[data-copy="datasheet"]').parentElement;
  row.insertAdjacentHTML('afterbegin', '<button class="water" id="dsSend">Send masses to post-lab</button>');
  $('#dsSend').addEventListener('click', () => {
    const map = { ba: ['in_mBA', 'in_wBA'], n: ['in_mN', 'in_wN'], c: ['in_mCC', 'in_wCC'] }; let sent = 0;
    DS_COLS.forEach(([c]) => { const m = num(`ds_mass_${c}`), wb = num(`ds_wb_${c}`); if (Number.isFinite(m)) { setVal(map[c][0], m); sent++; } if (Number.isFinite(wb)) { setVal(map[c][1], wb); sent++; } });
    $('[data-copied="datasheet"]').textContent = sent ? `Sent ${sent} values to the post-lab step checks.` : 'Fill in the sample and wire masses first.';
  });
}

/* =====================================================================
   POST-LAB
   ===================================================================== */
let azData = null, azName = '';
function analyzer() {
  const src = $('#azSrc'); if (virtualRun) src.querySelector('option[value="virtual"]').disabled = false;
  const sPre = $('#azPre'), sPost = $('#azPost'), sFrac = $('#azFrac');
  let last = null;
  const load = (d, name, tag) => {
    azData = d; azName = name; $('#azTag').textContent = tag;
    const n = d.t.length; sPre.max = n - 5; sPost.max = n - 3;
    const b = autoBounds(d); sPre.value = b.iF; sPost.value = b.iP; render();
  };
  const render = () => {
    if (!azData) return;
    let iF = +sPre.value, iP = +sPost.value;
    if (iP <= iF + 1) { iP = iF + 2; sPost.value = iP; }
    const f = analyze(azData, iF, iP, +sFrac.value); last = f;
    $('#azPreV').textContent = (f.tf / 60).toFixed(2) + ' min'; $('#azPostV').textContent = (f.tp / 60).toFixed(2) + ' min'; $('#azFracV').textContent = Math.round(+sFrac.value * 100) + ' %';
    plot($('#azChart'), azData, { fit: f, label: `Thermogram: ${azName}` });
    $('#azStats').innerHTML = [
      ['Fore drift', (f.pre.b * 60 * 1000).toFixed(2) + ' mK/min'], ['After drift', (f.post.b * 60 * 1000).toFixed(2) + ' mK/min'],
      ['Extrapolation time', (f.tb / 60).toFixed(2) + ' min'], ['Uncorrected rise', f.raw.toFixed(4) + ' K'], ['Adiabatic ΔT', f.dT.toFixed(4) + ' K']
    ].map(([a, b]) => `<div class="stat"><span>${a}</span><b>${b}</b></div>`).join('');
  };
  [sPre, sPost, sFrac].forEach(s => s.addEventListener('input', render));
  src.addEventListener('change', () => {
    const v = src.value; $('#azPasteWrap').hidden = v !== 'paste'; $('#azUnitWrap').hidden = v !== 'paste'; fb($('#azFb'), '', '');
    if (v === 'sim-ba') load(exampleRun('ba'), 'benzoic acid (simulated)', 'simulated data');
    if (v === 'sim-n') load(exampleRun('n'), 'naphthalene (simulated)', 'simulated data');
    if (v === 'sim-c') load(exampleRun('c'), 'CCDT (simulated)', 'simulated data');
    if (v === 'virtual' && virtualRun) load(virtualRun, `virtual ${CMP[virtualRun.key].name} run`, 'your virtual run');
  });
  $('#azLoad').addEventListener('click', () => {
    const mult = $('#azUnit').value === 'min' ? 60 : 1, t = [], T = [];
    $('#azPaste').value.split(/\r?\n/).forEach(line => { const p = line.trim().split(/[\s,;]+/).map(Number); if (p.length >= 2 && Number.isFinite(p[0]) && Number.isFinite(p[1])) { t.push(p[0] * mult); T.push(p[1]); } });
    if (t.length < 20) return fb($('#azFb'), 'bad', `Found only ${t.length} valid rows. Paste at least 20 rows with time and temperature in two numeric columns.`);
    const t0 = t[0]; load({ t: t.map(x => x - t0), T }, 'your data', 'your data');
    fb($('#azFb'), 'ok', `Loaded ${t.length} points. Check that the orange fit lines cover only the straight parts of the curve, and adjust the sliders if they don't.`);
  });
  $$('[data-send]').forEach(b => b.addEventListener('click', () => {
    if (!last) return; const id = { ba: 'in_dTBA', n: 'in_dTN', c: 'in_dTCC' }[b.dataset.send];
    setVal(id, last.dT.toFixed(4)); fb($('#azFb'), 'ok', `ΔT = ${last.dT.toFixed(4)} K sent to ${b.textContent}.`);
  }));
  load(exampleRun('ba'), 'benzoic acid (simulated)', 'simulated data');
}
function refs() {
  const mBA = num('in_mBA'), wBA = num('in_wBA') || 0, dTBA = num('in_dTBA');
  const C = (mBA * BA_DU * 1000 + wBA * WIRE * 1000) / dTBA;
  const CnoW = (mBA * BA_DU * 1000) / dTBA;
  const U = (m, w, dT, M) => -(C * dT - w * WIRE * 1000) / m * M / 1000;
  const UnoW = (m, dT, M) => -(C * dT) / m * M / 1000;
  const mN = num('in_mN'), wN = num('in_wN') || 0, dTN = num('in_dTN');
  const mC = num('in_mCC'), wC = num('in_wCC') || 0, dTC = num('in_dTCC');
  const UN = U(mN, wN, dTN, 128.17), UC = U(mC, wC, dTC, 162.28);
  const RT = R * 298.15 / 1000;
  const HN = UN - 2 * RT, HC = UC - 4.5 * RT;
  const fN = 10 * num('in_fCO2') + 4 * num('in_fH2O') - HN;
  const vC = num('in_tr') * (num('in_Tb') + 273.15) / 1000;
  const cB = num('in_cB'), vB = num('in_vB'), cX = num('in_cX'), vX = num('in_vX');
  const Er = (cB - vB) + (cX - vX) - (HC - vC);
  return { C, CnoW, UN, UC, HN, HC, fN, vC, Er, cB, vB, cX, vX, RT, mN, dTN, mC, dTC, mBA, dTBA, UnoWN: UnoW(mN, dTN, 128.17), UnoWC: UnoW(mC, dTC, 162.28) };
}
const near = (a, b, rel = 0.004, abs = 0) => Number.isFinite(a) && Number.isFinite(b) && Math.abs(a - b) <= Math.max(Math.abs(b) * rel, abs);
function checkStep(key) {
  const r = refs(), out = $('#fb_' + key);
  const signMsg = ['C', 'UN', 'UCC', 'HN', 'HCC'].includes(key) ? 'The magnitude is right but the sign is wrong. Combustion releases heat, so ΔcU and ΔcH are negative, while C is positive.' : 'The magnitude is right but the sign is wrong. Recheck the order of subtraction and the sign of each term.';
  const generic = (a, ref, rel, abs, extra = []) => {
    if (!Number.isFinite(ref)) return fb(out, 'warn', 'Some inputs above are missing or invalid. Fill them in first.');
    if (!Number.isFinite(a)) return fb(out, 'warn', 'Enter your answer first.');
    if (near(a, ref, rel, abs)) return fb(out, 'ok', '✓ Consistent with your inputs.');
    for (const [cond, msg] of extra) if (cond) return fb(out, 'warn', msg);
    if (near(a, -ref, rel, abs)) return fb(out, 'warn', signMsg);
    if (near(a, ref * 1000, rel) || near(a, ref / 1000, rel)) return fb(out, 'warn', 'Off by a factor of 1000. Check your J/kJ or g/mg conversions.');
    const off = Math.abs((a - ref) / ref) * 100;
    fb(out, 'bad', `Not consistent with your inputs (off by about ${off < 1 ? off.toFixed(2) : off.toFixed(1)} %). Recheck the arithmetic and units step by step.`);
  };
  const a = id => num(id);
  switch (key) {
    case 'C': return generic(a('ans_C'), r.C, 0.004, 0, [[near(a('ans_C'), r.CnoW, 0.002), 'Close, but this looks like C without the fuse-wire term. Add m<sub>wire</sub>·5.857 kJ g⁻¹ to the heat released.']]);
    case 'UN': return generic(a('ans_UN'), r.UN, 0.004, 0, [[near(a('ans_UN'), r.UnoWN, 0.002), 'You didn\'t subtract the heat from the fuse wire before dividing by the moles of sample.'], [near(a('ans_UN'), r.UN / 128.17, 0.01), 'That looks like a per-gram value. Multiply by the molar mass to get kJ mol⁻¹.']]);
    case 'UCC': return generic(a('ans_UCC'), r.UC, 0.004, 0, [[near(a('ans_UCC'), r.UnoWC, 0.002), 'You didn\'t subtract the heat from the fuse wire before dividing by the moles of sample.'], [near(a('ans_UCC'), r.UC / 162.28, 0.01), 'That looks like a per-gram value. Multiply by the molar mass to get kJ mol⁻¹.']]);
    case 'HN': case 'HCC': {
      const dnId = key === 'HN' ? 'ans_dnN' : 'ans_dnCC', dn = key === 'HN' ? -2 : -4.5, U = key === 'HN' ? r.UN : r.UC, H = key === 'HN' ? r.HN : r.HC, ans = a(key === 'HN' ? 'ans_HN' : 'ans_HCC');
      const sd = a(dnId);
      if (!Number.isFinite(sd)) return fb(out, 'warn', 'Enter Δn<sub>gas</sub> as well.');
      if (Math.abs(sd - dn) > 1e-6) return fb(out, 'warn', 'Check Δn<sub>gas</sub>. Count only gases: moles of CO₂ produced minus moles of O₂ consumed, per mole of fuel. H₂O(l) doesn\'t count.');
      return generic(ans, H, 0.0004, 0.6, [[near(ans, U, 0.0002, 0.3), 'This equals ΔcU. You still need to add Δn<sub>gas</sub>RT.'], [near(ans, U - dn * r.RT, 0.0002, 0.3), 'Check the sign of the correction: ΔH = ΔU + Δn<sub>gas</sub>RT, and with Δn<sub>gas</sub> &lt; 0, ΔH is more negative than ΔU.'], [near(ans, U + dn * R * 298.15, 0.001), 'R·T is in J mol⁻¹. Convert it to kJ before adding it to ΔU.']]);
    }
    case 'fN': return generic(a('ans_fN'), r.fN, 0, 1.5, [[near(a('ans_fN'), 10 * num('in_fCO2') + 4 * num('in_fH2O') + r.HN, 0, 1.5), 'You added ΔcH instead of subtracting it. ΔcH = ΣΔfH(products) − ΔfH(naphthalene).'], [near(a('ans_fN'), num('in_fCO2') + num('in_fH2O') - r.HN, 0, 1.5), 'Multiply the formation enthalpies by their stoichiometric coefficients (10 CO₂, 4 H₂O).']]);
    case 'vCC': return generic(a('ans_vCC'), r.vC, 0.004, 0, [[near(a('ans_vCC'), num('in_tr') * num('in_Tb') / 1000, 0.004), 'Use the boiling point in kelvin, not °C.']]);
    case 'R': {
      if (![r.cB, r.vB, r.cX, r.vX].every(Number.isFinite)) return fb(out, 'warn', 'Enter the literature values for benzene and cyclohexane first, and cite your source in the report.');
      const ans = a('ans_R');
      const res = generic(ans, r.Er, 0, 2, [
        [near(ans, r.cB + r.cX - r.HC, 0, 2), 'You left out the vaporization corrections. Convert every liquid ΔcH to the gas phase first.'],
        [near(ans, (r.cB + r.vB) + (r.cX + r.vX) - (r.HC + r.vC), 0, 2), 'Check the sign of the vaporization terms: ΔcH(g) = ΔcH(l) − ΔvapH.'],
        [near(ans, (r.cB - r.vB) + (r.cX - r.vX) - r.HC, 0, 2), 'You didn\'t correct CCDT to the gas phase. Subtract ΔvapH(CCDT) from its ΔcH(l) too.']]);
      if (near(ans, r.Er, 0, 2)) fb(out, 'ok', `✓ Consistent with your inputs. Now compare it with ≈ 150 kJ mol⁻¹ from hydrogenation data, and use the uncertainty explorer to judge whether the difference is significant.`);
      return res;
    }
  }
}
function uncertainty() {
  const ids = ['uM', 'uT', 'uV', 'uL'];
  const upd = () => {
    const sm = +$('#uM').value / 1000, sT = +$('#uT').value, sV = +$('#uV').value / 100, sL = +$('#uL').value;
    $('#uMv').textContent = (+$('#uM').value).toFixed(2); $('#uTv').textContent = sT.toFixed(3); $('#uVv').textContent = (sV * 100).toFixed(2); $('#uLv').textContent = sL.toFixed(1);
    const r = refs();
    const mBA = r.mBA > 0 ? r.mBA : 0.18, dTBA = r.dTBA > 0 ? r.dTBA : 2.13, mC = r.mC > 0 ? r.mC : 0.10, dTC = r.dTC > 0 ? r.dTC : 2.03;
    const H = Number.isFinite(r.HC) ? Math.abs(r.HC) : 7330;
    const relCal = Math.hypot(sm / mBA, sT / dTBA);
    const parts = [
      ['Sample mass', (sm / mC) * H], ['ΔT, CCDT run', (sT / dTC) * H], ['Calibration C', relCal * H],
      ['Bath volume', Math.SQRT2 * sV * H], ['Literature values', Math.SQRT2 * sL], ['Trouton estimate', 3 * ((num('in_Tb') || 231) + 273.15) / 1000]
    ];
    const sH = Math.hypot(parts[0][1], parts[1][1], parts[2][1], parts[3][1]);
    const sE = Math.hypot(...parts.map(p => p[1]));
    const E = Number.isFinite(r.Er) ? r.Er : 150;
    const mx = Math.max(...parts.map(p => p[1]));
    $('#uncBars').innerHTML = '<div class="bar" style="font-family:var(--mono);font-size:.75rem;color:var(--muted)"><span>contribution to σ(E<sub>res</sub>)</span><span></span><span style="text-align:right">kJ mol⁻¹</span></div>' +
      parts.map(([n, v], i) => `<div class="bar"><span>${n}</span><i class="${i >= 4 ? 'em' : ''}" style="width:${(v / mx * 100).toFixed(1)}%"></i><b>${v.toFixed(1)}</b></div>`).join('');
    $('#uncStats').innerHTML = [
      ['σ(ΔcH, CCDT)', `${sH.toFixed(1)} kJ mol⁻¹`], ['Relative σ(ΔcH)', `${(sH / H * 100).toFixed(2)} %`],
      ['σ(E_res)', `${sE.toFixed(1)} kJ mol⁻¹`], [`Relative σ(E_res)${Number.isFinite(r.Er) ? '' : ', using 150'}`, `${(sE / Math.abs(E) * 100).toFixed(0)} %`]
    ].map(([a, b]) => `<div class="stat"><span>${a}</span><b>${b}</b></div>`).join('');
  };
  ids.forEach(id => $('#' + id).addEventListener('input', upd));
  document.addEventListener('input', e => { if (e.target.id && e.target.id.startsWith('in_')) upd(); });
  upd();
}
const RUBRIC = [
  ['Part 1 · Thermograms (24)', ['Temperature vs time plots for benzoic acid, naphthalene and CCDT (6)', 'Axis labels with units on both axes (4)', 'Figure captions (4)', 'Linear fits to the fore and after periods shown for benzoic acid (5)', 'ΔT reported for benzoic acid, naphthalene and CCDT (5)']],
  ['Part 2 · Calorimeter heat capacity (9)', ['Δcu for benzoic acid and the fuse wire stated with units (2)', 'Total energy: E = m<sub>BA</sub>Δcu<sub>BA</sub> + m<sub>wire</sub>Δcu<sub>wire</sub> (3)', 'C = −E/ΔT, with a value around 2000–2500 J K⁻¹ (4)']],
  ['Part 3 · Heats of combustion (18)', ['Fuse-wire correction applied (2)', 'E = −C·ΔT (3)', 'ΔcU = E/n, with the molar mass of naphthalene (3)', 'ΔcH = ΔcU + Δn<sub>gas</sub>RT, with Δn<sub>gas</sub> shown (4)', 'Same analysis for CCDT, with M and Δn<sub>gas</sub> (2)', 'Uncertainty in ΔcH propagated from m, ΔT, C and bath volume (4)']],
  ['Part 4 · ΔfH° of naphthalene (19)', ['Balanced combustion equation with phases (3)', 'Literature ΔfH° of CO₂(g) and H₂O(l), cited (4)', 'Hess\'s law expression for ΔfH°(C₁₀H₈) (5) and value (1)', 'Uncertainty in ΔfH° (3)', 'Comparison with the literature value, with explanation (3)']],
  ['Part 5 · Resonance energy (30)', ['Trouton\'s rule estimate of ΔvapH(CCDT), using T<sub>b</sub> in K (4)', 'Literature ΔvapH for benzene and cyclohexane, cited (6)', 'E<sub>res</sub> expression with every component, phases corrected (10)', 'E<sub>res</sub> value (2)', 'Uncertainty in E<sub>res</sub> (3)', 'Comparison with literature (≈ 150 kJ mol⁻¹), the size of the uncertainty, a precision comment and an explanation (5)']]
];
function rubric() {
  let n = 0;
  $('#rubric').innerHTML = RUBRIC.map(([h, items]) => `<div class="phase"><h3 style="font-size:1.02rem">${h}</h3><ul class="steps">${items.map(s => { const id = 'rub_' + (n++); return `<li><input type="checkbox" id="${id}" data-save><label for="${id}">${s}</label></li>`; }).join('')}</ul></div>`).join('');
  $('#rubric').addEventListener('change', e => { const li = e.target.closest('li'); if (li) li.classList.toggle('done', e.target.checked); });
  rubric.upd = () => $$('#rubric input').forEach(b => b.closest('li').classList.toggle('done', b.checked));
}
const DISC = [
  ['1', 7, 'How much work (in joules) could be done under isobaric conditions by the combustion of one mole of (a) benzoic acid and (b) naphthalene? Assume an irreversible process.', 'At constant external pressure, w = −P<sub>ext</sub>ΔV. Treat the gases as ideal and count only gas moles, so PΔV = Δn<sub>gas</sub>RT. State your temperature, and be clear about the sign: is work done <i>by</i> the system or <i>on</i> it?'],
  ['2', 10, 'Compare the bomb calorimeter with the combustion chamber in a car engine. How is the energy released by combustion used differently in the two cases?', 'Rigid walls versus a moving piston. Where does the energy end up in each: as heat, as PV work, or lost to the surroundings?'],
  ['3', 16, 'Propose an experiment to find the resonance stabilization energy of naphthalene. Which compounds would you need to burn for an analysis like the one for benzene? Write balanced combustion equations and explain how these measurements would give a combustion energy to compare with naphthalene\'s, isolating the resonance energy. How should naphthalene\'s resonance energy compare with benzene\'s?', 'By analogy with C₁₂H₁₈ = C₆H₆ + C₆H₁₂, look for partly or fully hydrogenated bicyclic C₁₀ compounds whose combined bonds match a non-resonant "naphthalene". Revisit pre-lab question 5. Is the resonance energy per ring the same?'],
  ['4', 7, 'Suppose there is not enough oxygen in the bomb. How would the chemistry of combustion change? Write a suitable chemical equation. For naphthalene, how would this affect the energy released per gram of reactant?', 'Which carbon oxide forms when oxygen is limited? Compare its ΔfH° with that of CO₂. Why is this a health hazard in a home furnace?'],
  ['5', 10, 'Suppose you used an equal mass of linseed oil instead of a water bath. Which property of the bath liquid determines the temperature rise you measure? Predict the temperature rise for your benzoic acid standard in a linseed-oil bath. Which liquid would let you determine the enthalpy of combustion more accurately, and why?', 'q = m·c<sub>s</sub>·ΔT. Look up the specific heat capacity of linseed oil and cite it. Think about how a larger ΔT affects relative precision, and about any practical problems with oil (stirring, viscosity, heat transfer).']
];
function discQs() {
  $('#discQs').innerHTML = `<div class="panel"><div class="grid2"><div class="field"><label for="po_name">Name</label><input type="text" id="po_name" data-save autocomplete="name"></div><div class="field"><label for="po_sec">Lab section / date</label><input type="text" id="po_sec" data-save></div></div></div>` +
    DISC.map(([n, p, q, h]) => `<div class="q"><div class="q-head"><h4>Question ${n}</h4><span class="pts">${p} pts</span></div><p>${q}</p><label for="po_q${n}" class="sr">Answer to discussion question ${n}</label><textarea id="po_q${n}" data-save placeholder="Your answer…"></textarea><details class="hint"><summary>Hint</summary>${h}</details></div>`).join('');
}

/* ---------- copy answers ---------- */
const strip = h => { const d = document.createElement('div'); d.innerHTML = h; return d.textContent; };
function buildCopy(which) {
  const v = id => (document.getElementById(id)?.value || '').trim();
  if (which === 'prelab') {
    return `BOMB CALORIMETRY · PRE-LAB\nName: ${v('pl_name')}\nSection/date: ${v('pl_sec')}\n\n` + PRELAB_QS.map(([n, p, q]) => `Q${n} (${p} pts). ${strip(q)}\n${v('pl_q' + n) || '[no answer]'}\n`).join('\n');
  }
  if (which === 'postlab') {
    const rows = [['C (J/K)', 'ans_C'], ['ΔcU naphthalene (kJ/mol)', 'ans_UN'], ['Δn gas naphthalene', 'ans_dnN'], ['ΔcH naphthalene (kJ/mol)', 'ans_HN'], ['ΔcU CCDT (kJ/mol)', 'ans_UCC'], ['Δn gas CCDT', 'ans_dnCC'], ['ΔcH CCDT (kJ/mol)', 'ans_HCC'], ['ΔfH naphthalene (kJ/mol)', 'ans_fN'], ['ΔvapH CCDT (kJ/mol)', 'ans_vCC'], ['E_res (kJ/mol)', 'ans_R']];
    const inputs = [['m BA (g)', 'in_mBA'], ['wire BA (g)', 'in_wBA'], ['ΔT BA (K)', 'in_dTBA'], ['m naph (g)', 'in_mN'], ['wire naph (g)', 'in_wN'], ['ΔT naph (K)', 'in_dTN'], ['m CCDT (g)', 'in_mCC'], ['wire CCDT (g)', 'in_wCC'], ['ΔT CCDT (K)', 'in_dTCC']];
    return `BOMB CALORIMETRY · POST-LAB\nName: ${v('po_name')}\nSection/date: ${v('po_sec')}\n\nINPUT DATA\n` + inputs.map(([l, id]) => `  ${l}: ${v(id)}`).join('\n') + `\n\nRESULTS\n` + rows.map(([l, id]) => `  ${l}: ${v(id) || '—'}`).join('\n') + `\n\nDISCUSSION\n` + DISC.map(([n, p, q]) => `Q${n} (${p} pts). ${strip(q)}\n${v('po_q' + n) || '[no answer]'}\n`).join('\n');
  }
  if (which === 'datasheet') {
    return 'BOMB CALORIMETRY · DATA SHEET\n' + ['Quantity', ...DS_COLS.map(c => c[1])].join('\t') + '\n' + DS_ROWS.map(([k, l]) => [l, ...DS_COLS.map(([c]) => v(`ds_${k}_${c}`))].join('\t')).join('\n');
  }
  return '';
}
document.addEventListener('click', e => {
  const b = e.target.closest('[data-copy]'); if (!b) return;
  const which = b.dataset.copy, text = buildCopy(which), note = $(`[data-copied="${which}"]`);
  const fallback = () => { const ta = $('#copyFallback-' + which); ta.classList.remove('sr'); ta.removeAttribute('aria-hidden'); ta.value = text; ta.style.marginTop = '10px'; ta.focus(); ta.select(); note.textContent = 'Your browser blocked automatic copying. The text is selected below. Press Ctrl+C (⌘C on a Mac) to copy it.'; };
  try { navigator.clipboard.writeText(text).then(() => { note.textContent = 'Copied. Paste it into your course submission.'; }, fallback); } catch (err) { fallback(); }
});

/* =====================================================================
   INIT
   ===================================================================== */
function init() {
  initCal();
  molecules(); balancer(); dnWidget(); calibWidget(); hess(); sorter(); quiz(); prelabQs();
  chargeCheck(); thermoFig(); simulator(); procedure(); dataSheet();
  analyzer(); rubric(); discQs();
  bindSaves();
  procedure.upd(); dataSheet.upd(); rubric.upd();
  setTimeout(() => $('#balancer').dispatchEvent(new Event('input')), 0);
  $$('[data-check]').forEach(b => b.addEventListener('click', () => checkStep(b.dataset.check)));
  uncertainty();
  route();
}
init();
})();
