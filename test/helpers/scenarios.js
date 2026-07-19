// Scenario driving + figure snapshotting for the golden suite.
//
// A "scenario" is a set of slider positions + y/x axis choices (+ optional ?ext
// decade shift) applied to a page, after which update() runs and we snapshot the
// resulting figure data (curves, pie, axes, data table). These snapshots are the
// golden reference that every later refactor step must reproduce byte-for-byte.

import { loadSim } from './loadSim.js';

// Radio conventions (from the pages):
//   ligand/homodimer/ligands y-axis: 1=absolute 2=relative ; x-axis: 3 / 4
//   receptors y-axis: 5=abs-log 6=abs-lin 7=specificity   ; x-axis: 3 / 4
export const SCENARIOS = {
  ligand: [
    { name: 'abs-total',   yradio: 1, xradio: 3 },
    { name: 'rel-total',   yradio: 2, xradio: 3 },
    { name: 'abs-free',    yradio: 1, xradio: 4 },
    { name: 'rel-free',    yradio: 2, xradio: 4 },
    { name: 'alt-sliders', yradio: 1, xradio: 3, sliders: { 5: 180, 3: 60, 7: 150 } },
    { name: 'ext-decade',  yradio: 1, xradio: 3, ext: true, ds: 2 },
  ],
  homodimer: [
    { name: 'abs-total',   yradio: 1, xradio: 3 },
    { name: 'rel-total',   yradio: 2, xradio: 3 },
    { name: 'abs-free',    yradio: 1, xradio: 4 },
    { name: 'rel-free',    yradio: 2, xradio: 4 },
    { name: 'alt-sliders', yradio: 1, xradio: 3, sliders: { 5: 200, 7: 120 } },
  ],
  ligands: [
    { name: 'abs-ligandX', yradio: 1, xradio: 4 },
    { name: 'rel-ligandX', yradio: 2, xradio: 4 },
    { name: 'abs-proteinX', yradio: 1, xradio: 3 },
    { name: 'rel-proteinX', yradio: 2, xradio: 3 },
    { name: 'alt-sliders', yradio: 1, xradio: 4, sliders: { 5: 150, 10: 100, 3: 120, 7: 120, 9: 150 } },
  ],
  receptors: [
    { name: 'absLog-total', yradio: 5, xradio: 3 },
    { name: 'absLin-total', yradio: 6, xradio: 3 },
    { name: 'spec-total',   yradio: 7, xradio: 3 },
    { name: 'absLog-free',  yradio: 5, xradio: 4 },
    { name: 'spec-free',    yradio: 7, xradio: 4 },
    { name: 'alt-sliders',  yradio: 5, xradio: 3, sliders: { 5: 150, 10: 100, 3: 120, 7: 120, 9: 150 } },
  ],
};

/** Load a page, apply a scenario, run update(), and return the live window. */
export function runScenario(pageKey, sc) {
  const w = loadSim(pageKey, { ext: !!sc.ext });

  if (sc.ext && sc.ds != null) {
    w.document.getElementById('ext_dsinput').value = String(sc.ds);
    w.decadeshift_changed(true);   // set decadeshift, defer update
  }

  if (sc.sliders) {
    for (const [idx, val] of Object.entries(sc.sliders)) {
      w.document.getElementById('slider' + idx).value = String(val);
      w.slider_input(Number(idx), true);   // recompute global, no update
    }
  }

  if (sc.yradio) w.radio_input(sc.yradio, true);
  if (sc.xradio) w.radio_input(sc.xradio, true);

  w.update();
  return w;
}

function readDataTable(w) {
  const labels = w.datalabels[w.appmode] || [];
  const rows = [];
  for (let i = 0; i < labels.length; i++) {
    const cell = id => { const e = w.document.getElementById('data' + i + id); return e ? e.innerHTML : null; };
    rows.push({ a: cell('a'), b: cell('b'), c: cell('c') });
  }
  return rows;
}

/** Capture the figure-defining state after update(). Order of keys is stable. */
export function snapshot(w) {
  return {
    globals: {
      S_0: w.S_0, E_0: w.E_0, K_D: w.K_D, K_D2: w.K_D2, Q_0: w.Q_0,
      scale_absolute: w.scale_absolute, xscale_alternative: w.xscale_alternative,
      decadeshift: w.decadeshift,
      xmin: w.xmin, xmax: w.xmax, ymin: w.ymin, ymax: w.ymax,
      xaxistype: w.xaxistype, yaxistype: w.yaxistype,
    },
    curves: w.curves.map(c => Array.isArray(c) ? c.map(p => [p.x, p.y]) : null),
    pie: w.pie.map(p => p.value),
    labels: w.labels.map(l => (l == null ? null : l)),
    colours: w.colours.slice(),
    legends: w.legends.map(l => (l == null ? null : l)),
    dataTable: readDataTable(w),
  };
}

const nonFinite = (k, v) => (typeof v === 'number' && !Number.isFinite(v) ? '@@' + String(v) : v);

/** Deterministic, human-readable JSON; non-finite numbers encoded for exact equality. */
export function stableStringify(obj) {
  return JSON.stringify(obj, nonFinite, 2);
}

/** Same encoding, no indentation — for large goldens (curves) where size matters. */
export function compactStringify(obj) {
  return JSON.stringify(obj, nonFinite);
}

/** Build the full curves-golden object: { model: { scenario: snapshot } }. */
export function buildCurvesGolden() {
  const out = {};
  for (const [pageKey, list] of Object.entries(SCENARIOS)) {
    out[pageKey] = {};
    for (const sc of list) {
      out[pageKey][sc.name] = snapshot(runScenario(pageKey, sc));
    }
  }
  return out;
}
