// Tolerant structural comparison for the golden gate.
//
// Golden JSON is captured on one machine but the suite runs on several (local
// macOS + Linux CI), so numbers are compared within a relative tolerance rather
// than bit-exactly. Two effects require this:
//   1. Transcendental functions (Math.pow/log10/log) differ by ~1 ULP (~1e-16 rel)
//      between platforms' math libraries — negligible.
//   2. A handful of curve points sit in a catastrophic-cancellation regime: the
//      free-protein species [P] = E_0 - v, where v -> E_0 at high ligand, loses
//      most of its significant figures (see the numerical-stability note in
//      calculate_ligand_total / REFACTOR_PLAN.md). Those tail points are near the
//      noise floor and are NOT bit-reproducible across platforms; observed
//      cross-platform delta is up to ~1.7e-4 relative on ~11 isolated points.
//
// REL = 1e-3 absorbs (2) with ~6x margin while still catching any real behaviour
// change: a correct refactor is bit-identical on a given platform (~1e-15 at
// worst from re-association), and a genuine bug (wrong formula/constant/branch)
// shifts values by whole percent — orders of magnitude above this tolerance.
// The analytic invariants (invariants.test.js) remain tight (1e-12) as a second,
// well-conditioned line of defence.
//
// Strings (labels, colours, formatted data-table cells) and booleans are compared
// exactly. Structure (keys, array lengths, null/undefined) must match exactly.

const REL = 1e-3;   // relative tolerance (see note above re: cancellation tail)
const ABS = 1e-12;  // absolute floor, for values near zero

function numClose(a, b) {
  if (Number.isNaN(a) && Number.isNaN(b)) return true;
  if (!Number.isFinite(a) || !Number.isFinite(b)) return a === b;   // ±Infinity must match
  return Math.abs(a - b) <= ABS + REL * Math.max(Math.abs(a), Math.abs(b));
}

/** Parse a golden file, decoding the non-finite sentinels written by the serializers. */
export function parseGolden(text) {
  return JSON.parse(text, (k, v) => {
    if (typeof v === 'string' && v.length > 2 && v[0] === '@' && v[1] === '@') {
      const s = v.slice(2);
      if (s === 'NaN') return NaN;
      if (s === 'Infinity') return Infinity;
      if (s === '-Infinity') return -Infinity;
    }
    return v;
  });
}

/**
 * Walk `actual` vs `expected`, collecting up to `limit` human-readable mismatch
 * descriptions (with JSON-ish paths). Empty result means they match within tolerance.
 */
export function findMismatches(actual, expected, path = '', out = [], limit = 25) {
  if (out.length >= limit) return out;

  if (typeof expected === 'number' || typeof actual === 'number') {
    if (typeof actual !== typeof expected || !numClose(actual, expected))
      out.push(`${path}: expected ${expected}, got ${actual}`);
    return out;
  }

  if (Array.isArray(expected) || Array.isArray(actual)) {
    if (!Array.isArray(expected) || !Array.isArray(actual) || actual.length !== expected.length) {
      out.push(`${path}: array shape differs (len ${actual && actual.length} vs ${expected && expected.length})`);
      return out;
    }
    for (let i = 0; i < expected.length && out.length < limit; i++)
      findMismatches(actual[i], expected[i], `${path}[${i}]`, out, limit);
    return out;
  }

  if (expected && typeof expected === 'object') {
    if (!actual || typeof actual !== 'object') { out.push(`${path}: expected object, got ${actual}`); return out; }
    for (const k of Object.keys(expected)) {
      if (out.length >= limit) break;
      if (!(k in actual)) { out.push(`${path}.${k}: missing`); continue; }
      findMismatches(actual[k], expected[k], `${path}.${k}`, out, limit);
    }
    for (const k of Object.keys(actual))
      if (!(k in expected)) out.push(`${path}.${k}: unexpected key`);
    return out;
  }

  // strings, booleans, null
  if (actual !== expected)
    out.push(`${path}: expected ${JSON.stringify(expected)}, got ${JSON.stringify(actual)}`);
  return out;
}
