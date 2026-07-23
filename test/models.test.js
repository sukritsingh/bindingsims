// Model-contract suite (refactor/03): generic invariants that every registered
// model must satisfy. New models inherit this coverage automatically by being
// exported from models/index.js — the guardrail against registry drift.

import { test, expect, describe } from 'vitest';
import { models, modelByAppmode, datalabels } from '../models/index.js';

const entries = Object.entries(models);

test('registry is non-empty and keyed by descriptor id', () => {
  expect(entries.length).toBeGreaterThan(0);
  for (const [key, m] of entries) expect(m.id).toBe(key);
});

test('legacy appmode ids are unique and round-trip through modelByAppmode', () => {
  const seen = new Set();
  for (const [, m] of entries) {
    expect(Number.isInteger(m.legacyAppmode)).toBe(true);
    expect(seen.has(m.legacyAppmode)).toBe(false);
    seen.add(m.legacyAppmode);
    expect(modelByAppmode(m.legacyAppmode)).toBe(m);
  }
});

test('datalabels export mirrors each descriptor by legacy appmode', () => {
  for (const [, m] of entries) {
    expect(Array.isArray(m.datalabels)).toBe(true);
    // may be empty for rate models with no species table (e.g. inhibition)
    for (const label of m.datalabels) expect(typeof label).toBe('string');
    expect(datalabels[m.legacyAppmode]).toBe(m.datalabels);
  }
});

describe.each(entries)('descriptor: %s', (id, m) => {
  test('axisLabels returns {x,y} strings for every scale mode', () => {
    for (const xscale_alternative of [false, true]) {
      for (const scale_absolute of [0, 1, 2]) {
        const a = {
          xscale_alternative, scale_absolute, xmagnitude: 0, ymagnitude: -3,
          magnitude_string: (mag, no_trailing_space) =>
            !mag ? '' : '10^' + mag + (no_trailing_space === true ? '' : ' '),
        };
        const out = m.axisLabels(a);
        expect(typeof out.x).toBe('string');
        expect(typeof out.y).toBe('string');
      }
    }
  });

  test('fitSolve returns a {d,m} shape for both x-axis modes', () => {
    // m indexed by slider id (see expparams); representative micromolar values.
    const params = []; for (const idx of [3, 5, 7, 9, 10]) params[idx] = 1e-6;
    // The competing solvers read the ambient appmode globals (as in the app,
    // where they resolve off window); supply them for the isolated unit test.
    globalThis.appmode = m.legacyAppmode;
    globalThis.appmode_ligands = 0;
    try {
      for (const xscale_alternative of [false, true]) {
        const r = m.fitSolve(params, 1e-6, xscale_alternative);
        expect(Array.isArray(r.d)).toBe(true);
        expect(Array.isArray(r.m)).toBe(true);
        expect(r.d.length).toBe(r.m.length);
      }
    } finally {
      delete globalThis.appmode;
      delete globalThis.appmode_ligands;
    }
  });
});
