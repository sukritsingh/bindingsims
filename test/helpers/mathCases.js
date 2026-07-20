// Golden inputs for the pure math (G1) and the fitter (G3).
//
// G1 pins the arithmetic of every solver and the slider->value mapping (expval)
// against fixed input grids. G3 pins calculate_lsf: from noise-free synthetic data
// it records the recovered slider positions. Both are captured from the current,
// unmodified code and frozen.

import { loadSim } from './loadSim.js';

const CONC = [1e-9, 1e-7, 1e-6, 1e-5, 1e-3, 1e-1];

function grid3(fn) {
  const rows = [];
  for (const a of [1e-6, 1e-5])
    for (const b of [1e-8, 1e-6, 1e-4])
      for (const c of [1e-7, 1e-6, 1e-5])
        rows.push({ args: [a, b, c], out: fn(a, b, c) });
  return rows;
}

export function buildMathGolden() {
  const w = loadSim('ligands');           // has every solver; appmode = ligands
  const res = { r: v => (v && v.d ? { d: v.d, m: v.m } : v) };

  // expval — the slider contract. Production path (newdecadescale = true).
  const expTrue = {};
  for (const ds of [-2, 0, 3]) {
    w.decadeshift = ds;
    expTrue['ds' + ds] = [];
    for (let v = 0; v <= 240; v++) expTrue['ds' + ds].push(w.expval(v, 10, -9, -1));
  }
  w.decadeshift = 0;
  // Dev path (newdecadescale = false) — sampled.
  w.newdecadescale = false;
  const expFalse = {};
  for (const v of [0, 30, 60, 120, 180, 240]) expFalse['v' + v] = w.expval(v, 10, -9, -1);
  w.newdecadescale = true;

  const golden = {
    expval: { newdecadescale_true: expTrue, newdecadescale_false: expFalse },

    ligand_free:  grid3((E, S, K) => res.r(w.calculate_ligand_free(E, S, K))),
    ligand_total: grid3((E, S, K) => res.r(w.calculate_ligand_total(E, S, K))),

    homodimer:      CONC.flatMap(P => [1e-7, 1e-6, 1e-5].map(K => ({ args: [P, K], out: res.r(w.calculate_homodimer(P, K)) }))),
    homodimer_free: CONC.flatMap(P => [1e-7, 1e-6, 1e-5].map(K => ({ args: [P, K], out: res.r(w.calculate_homodimer_free(P, K)) }))),

    // 1 protein + 2 competing ligands (cubic). ligands-shaped output.
    ligands: (() => {
      const rows = [];
      for (const cA of [1e-6, 1e-5])
        for (const cB of [1e-7, 1e-5])
          for (const cC of [1e-7, 1e-5])
            for (const KD of [1e-6]) for (const KD2 of [1e-7, 1e-5])
              rows.push({ args: [cA, cB, cC, KD, KD2], out: res.r(w.calculate_ligands(cA, cB, cC, KD, KD2)) });
      return rows;
    })(),
  };

  // calculate_ligands_free depends on appmode for its output shape; capture both.
  const freeInputs = [];
  for (const A of [1e-7, 1e-6, 1e-5])
    for (const cB of [1e-6]) for (const cC of [1e-6])
      for (const KD of [1e-6]) for (const KD2 of [1e-6, 1e-5])
        freeInputs.push([A, cB, cC, KD, KD2]);

  w.appmode = w.appmode_ligands;
  golden.ligands_free_ligands = freeInputs.map(a => ({ args: a, out: res.r(w.calculate_ligands_free(...a)) }));
  w.appmode = w.appmode_receptors;
  golden.ligands_free_receptors = freeInputs.map(a => ({ args: a, out: res.r(w.calculate_ligands_free(...a)) }));
  w.appmode = w.appmode_ligands;

  // Cubic solvers on a few coefficient sets drawn from the competition polynomial.
  const cubicSets = [
    [1, 2e-5, 1e-11, -1e-17, 1e-6],
    [1, -3e-6, 2e-12, -1e-19, 5e-7],
    [1, 1e-4, 1e-9, -1e-15, 1e-5],
  ];
  golden.cubic = cubicSets.map(([q1, q2, q3, q4, a0]) => ({
    args: [q1, q2, q3, q4, a0],
    newton: w.solve_cubic_newton(q1, q2, q3, q4, a0),
    bisection: w.solve_cubic_bisection(q1, q2, q3, q4, 0, a0),
  }));

  return golden;
}

// ---- Fitter (G3) --------------------------------------------------------------

function setRadioGroup(w, name, id) {
  for (const el of w.document.getElementsByName(name)) el.checked = false;
  const chosen = w.document.getElementById(id);
  if (chosen) chosen.checked = true;
}

function runFit(cfg) {
  const w = loadSim(cfg.page);
  w.radio_input(cfg.yradio, true);
  w.radio_input(cfg.xradio, true);

  const data = cfg.data(w);
  w.document.getElementById('databox').value = data.map(d => d.x + ' ' + d.y).join('\n');
  w.data_changed(false);

  for (const idx of cfg.free) {
    const cb = w.document.getElementById('fixval' + idx);
    cb.disabled = false; cb.checked = true;
  }
  setRadioGroup(w, 'calcmode', 'calcmode' + cfg.calcmode);
  setRadioGroup(w, 'calcoption', 'calcoption' + cfg.calcoption);

  w.calculate_lsf();

  const sliders = {};
  for (const idx of cfg.free) sliders['slider' + idx] = w.document.getElementById('slider' + idx).value;
  return {
    sliders,
    status: w.document.getElementById('calcstatus').innerHTML,
    globals: { K_D: w.K_D, K_D2: w.K_D2, E_0: w.E_0, S_0: w.S_0, Q_0: w.Q_0 },
  };
}

// Synthetic data generators: noise-free, relative [PL]-type signal, true K_D = 1e-6.
const ligandData = w => [1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4]
  .map(x => { const d = w.calculate_ligand_total(w.E_0, x, 1e-6).d; return { x, y: d[0] / w.E_0 }; });

const homodimerData = w => [1e-8, 1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 1e-4]
  .map(x => { const d = w.calculate_homodimer(x, 1e-6).d; return { x, y: 2 * d[0] / (2 * d[0] + d[1]) }; });

const ligandsData = w => [1e-8, 1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 1e-4]
  .map(x => { const d = w.calculate_ligands(x, w.Q_0, w.S_0, 1e-6, w.K_D2).d; return { x, y: d[0] / (d[0] + d[1] + d[2]) }; });

export const FIT_CASES = {
  'ligand/calcmode1':   { page: 'ligand',    yradio: 2, xradio: 3, free: [7], calcmode: 1, calcoption: 0, data: ligandData },
  'ligand/calcmode2':   { page: 'ligand',    yradio: 2, xradio: 3, free: [7], calcmode: 2, calcoption: 0, data: ligandData },
  'ligand/calcmode3':   { page: 'ligand',    yradio: 2, xradio: 3, free: [7], calcmode: 3, calcoption: 0, data: ligandData },
  'homodimer/calcmode1': { page: 'homodimer', yradio: 2, xradio: 3, free: [7], calcmode: 1, calcoption: 0, data: homodimerData },
  'ligands/calcmode1':  { page: 'ligands',   yradio: 2, xradio: 3, free: [7], calcmode: 1, calcoption: 0, data: ligandsData },
  // calcmode 4 (extmode "precise iterative" — Jacobian + Gauss-Jordan covariance).
  // Added before extracting core/fit.js so the covariance path is golden-protected;
  // the 2-param case exercises the 2x2 matrix inversion. Recovered params land in
  // the globals (calcmode 4 writes ext_value + globals, not the slider positions).
  'ligand/calcmode4':    { page: 'ligand',    yradio: 2, xradio: 3, free: [7],    calcmode: 4, calcoption: 0, data: ligandData },
  'ligand/calcmode4-2p': { page: 'ligand',    yradio: 2, xradio: 3, free: [5, 7], calcmode: 4, calcoption: 0, data: ligandData },
};

export function buildFitGolden() {
  const out = {};
  for (const [name, cfg] of Object.entries(FIT_CASES)) out[name] = runFit(cfg);
  return out;
}
