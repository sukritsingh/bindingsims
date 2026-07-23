// Science tests for the enzyme-inhibition model (feature/07). These pin the
// corrected α/α′ forms and Cheng–Prusoff independently of the UI.

import { test, expect, describe } from 'vitest';
import { enzyme_rate, ic50, ic50_to_Ki, Ki_to_ic50 } from '../models/inhibition.js';

const base = { kcat: 100, E0: 1e-8, Km: 1e-5, Ki: 1e-6, KiP: 1e-6 };
const Vmax = base.kcat * base.E0;                 // 1e-6 mol l⁻¹ s⁻¹
const MECHS = ['none', 'competitive', 'uncompetitive', 'noncompetitive', 'mixed'];
const rel = (a, b) => Math.abs(a - b) <= 1e-12 * Math.max(1, Math.abs(a), Math.abs(b));

describe('reduces to plain Michaelis–Menten with no inhibitor', () => {
  test.each(MECHS)('mechanism %s at [I]=0', (mechanism) => {
    const p = { ...base, mechanism };
    for (const S of [1e-7, 1e-6, base.Km, 1e-4, 1e-2]) {
      const mm = Vmax * S / (base.Km + S);        // uninhibited MM
      expect(rel(enzyme_rate(p, S, 0), mm)).toBe(true);
    }
    // half-maximal at S = Km
    expect(rel(enzyme_rate(p, base.Km, 0), Vmax / 2)).toBe(true);
  });
});

test('competitive: Vmax unchanged at saturating S, apparent Km raised', () => {
  const p = { ...base, mechanism: 'competitive' };
  const I = 5e-6;                                  // α = 1 + I/Ki = 6
  expect(rel(enzyme_rate(p, 1e6, I), Vmax)).toBe(true);            // S → ∞ recovers Vmax
  // apparent Km = α·Km: half-maximal now at S = 6·Km
  expect(rel(enzyme_rate(p, 6 * base.Km, I), Vmax / 2)).toBe(true);
});

test('uncompetitive: Vmax and Km both lowered by α′', () => {
  const p = { ...base, mechanism: 'uncompetitive' };
  const I = 5e-6;                                  // α′ = 6
  expect(rel(enzyme_rate(p, 1e6, I), Vmax / 6)).toBe(true);        // apparent Vmax = Vmax/α′
  // apparent Km = Km/α′: half-maximal at S = Km/6
  expect(rel(enzyme_rate(p, base.Km / 6, I), Vmax / 6 / 2)).toBe(true);
});

test('non-competitive: Vmax lowered, apparent Km unchanged', () => {
  const p = { ...base, mechanism: 'noncompetitive' };
  const I = 5e-6;                                  // α = α′ = 6
  expect(rel(enzyme_rate(p, 1e6, I), Vmax / 6)).toBe(true);        // apparent Vmax = Vmax/6
  // still half-maximal at S = Km (Km unchanged)
  expect(rel(enzyme_rate(p, base.Km, I), Vmax / 6 / 2)).toBe(true);
});

test('non-competitive equals mixed with Ki′ = Ki', () => {
  const nc = { ...base, mechanism: 'noncompetitive' };
  const mx = { ...base, KiP: base.Ki, mechanism: 'mixed' };
  for (const S of [1e-7, 1e-5, 1e-3]) for (const I of [0, 1e-6, 1e-5])
    expect(rel(enzyme_rate(nc, S, I), enzyme_rate(mx, S, I))).toBe(true);
});

describe('Cheng–Prusoff: at [I]=IC50 the velocity is half the uninhibited value', () => {
  test.each(['competitive', 'uncompetitive', 'noncompetitive', 'mixed'])('%s', (mechanism) => {
    const p = { ...base, KiP: mechanism === 'mixed' ? 3e-6 : base.KiP, mechanism };
    for (const S of [3e-6, 1e-5, 5e-5]) {
      const half = enzyme_rate(p, S, 0) / 2;
      expect(rel(enzyme_rate(p, S, ic50(p, S)), half)).toBe(true);
    }
  });
});

test('classic competitive Cheng–Prusoff round-trips IC50 ↔ Ki', () => {
  const S = 2e-5, Km = base.Km, Ki = base.Ki;
  const back = ic50_to_Ki(Ki_to_ic50(Ki, S, Km), S, Km);
  expect(rel(back, Ki)).toBe(true);
});
