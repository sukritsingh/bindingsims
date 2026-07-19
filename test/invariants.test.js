// Analytic invariants — checks with a known-correct answer, independent of the
// golden snapshots. These encode the biophysics (mass conservation, half-saturation,
// parameter recovery) so a future change that breaks the science fails loudly.

import { test, expect, beforeAll } from 'vitest';
import { loadSim } from './helpers/loadSim.js';
import { FIT_CASES } from './helpers/mathCases.js';

let w;
beforeAll(() => { w = loadSim('ligands'); });   // exposes every solver; appmode = ligands

const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-30);

test('expval maps known slider positions exactly (production scale)', () => {
  w.decadeshift = 0;
  expect(w.expval(0, 10, -9, -1)).toBeCloseTo(1e-9, 20);
  expect(w.expval(30, 10, -9, -1)).toBeCloseTo(1e-8, 20);
  expect(w.expval(90, 10, -9, -1)).toBeCloseTo(1e-6, 20);
});

test('ligand (total form): mass is conserved for protein and ligand', () => {
  for (const E of [1e-6, 1e-5]) for (const S of [1e-8, 1e-6, 1e-4]) for (const K of [1e-7, 1e-6]) {
    const [PL, P, L] = w.calculate_ligand_total(E, S, K).d;
    expect(rel(PL + P, E)).toBeLessThan(1e-9);   // [PL] + [P] = E_0
    expect(rel(PL + L, S)).toBeLessThan(1e-9);   // [PL] + [L] = S_0
  }
});

test('ligand (free form): half saturation at [L] = K_D', () => {
  for (const E of [1e-6, 1e-5]) for (const K of [1e-7, 1e-6, 1e-5]) {
    const PL = w.calculate_ligand_free(E, K, K).d[0];
    expect(rel(PL, E / 2)).toBeLessThan(1e-12);
  }
});

test('homodimer: mass is conserved (2·[P2] + [P] = E_0)', () => {
  for (const E of [1e-6, 1e-5, 1e-4]) for (const K of [1e-7, 1e-6, 1e-5]) {
    const [P2, P] = w.calculate_homodimer(E, K).d;
    expect(rel(2 * P2 + P, E)).toBeLessThan(1e-9);
  }
});

test('competing ligands (cubic): mass is conserved for protein and both ligands', () => {
  for (const cA of [1e-6, 1e-5]) for (const cB of [1e-7, 1e-5]) for (const cC of [1e-7, 1e-5]) {
    const [AB, AC, A, B, C] = w.calculate_ligands(cA, cB, cC, 1e-6, 1e-5).d;
    expect(rel(AB + AC + A, cA)).toBeLessThan(1e-3);   // cubic solver tolerance ~1e-3
    expect(rel(AB + B, cB)).toBeLessThan(1e-3);
    expect(rel(AC + C, cC)).toBeLessThan(1e-3);
  }
});

test('cubic Newton and bisection agree on the competition root', () => {
  // Coefficients of the competition polynomial for a representative case.
  const cA = 1e-5, cB = 1e-6, cC = 1e-6, KD = 1e-6, KD2 = 1e-5;
  const t2 = KD2 + KD + cC + cB - cA;
  const t3 = KD2 * KD + cC * KD + cB * KD2 - cA * KD2 - cA * KD;
  const t4 = -cA * KD2 * KD;
  const n = w.solve_cubic_newton(1, t2, t3, t4, cA);
  const b = w.solve_cubic_bisection(1, t2, t3, t4, 0, cA);
  expect(rel(n, b)).toBeLessThan(1e-3);
  expect(n).toBeGreaterThanOrEqual(0);
  expect(n).toBeLessThanOrEqual(cA);
});

test('fitter recovers the true K_D (1e-6) from noise-free data', () => {
  // Independent re-run of the ligand two-pass fit case.
  const cfg = FIT_CASES['ligand/calcmode1'];
  const win = loadSim(cfg.page);
  win.radio_input(cfg.yradio, true);
  win.radio_input(cfg.xradio, true);
  const data = cfg.data(win);
  win.document.getElementById('databox').value = data.map(d => d.x + ' ' + d.y).join('\n');
  win.data_changed(false);
  const cb = win.document.getElementById('fixval7'); cb.disabled = false; cb.checked = true;
  win.document.getElementById('calcmode1').checked = true;
  win.document.getElementById('calcoption0').checked = true;
  win.calculate_lsf();
  expect(rel(win.K_D, 1e-6)).toBeLessThan(1e-6);
});
