// Golden gate: re-derive every snapshot from the CURRENT code and assert it equals
// the committed reference captured in step 00. This is what guarantees the refactor
// stays numerically identical. If a step legitimately changes behaviour, regenerate
// with `npm run golden` — never edit the JSON by hand.

import { test, expect } from 'vitest';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { buildCurvesGolden, stableStringify, compactStringify } from './helpers/scenarios.js';
import { buildMathGolden, buildFitGolden } from './helpers/mathCases.js';

const GOLDEN = path.resolve(path.dirname(fileURLToPath(import.meta.url)), 'golden');
const read = f => fs.readFileSync(path.join(GOLDEN, f), 'utf8').replace(/\n$/, '');

test('math golden (solvers + expval) unchanged', () => {
  expect(stableStringify(buildMathGolden())).toBe(read('math.golden.json'));
});

test('curves golden (figure data for all models/scenarios) unchanged', () => {
  expect(compactStringify(buildCurvesGolden())).toBe(read('curves.golden.json'));
});

test('fit golden (calculate_lsf recovery) unchanged', () => {
  expect(stableStringify(buildFitGolden())).toBe(read('fit.golden.json'));
});
