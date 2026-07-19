// Golden gate: re-derive every snapshot from the CURRENT code and assert it equals
// the committed reference captured in step 00, within a tight numeric tolerance
// (see helpers/compare.js — tolerance is for cross-platform float portability, not
// for masking behaviour change). This is what guarantees the refactor stays
// numerically identical. If a step legitimately changes behaviour, regenerate with
// `npm run golden` — never edit the JSON by hand.

import { test, expect } from 'vitest';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { buildCurvesGolden } from './helpers/scenarios.js';
import { buildMathGolden, buildFitGolden } from './helpers/mathCases.js';
import { parseGolden, findMismatches } from './helpers/compare.js';

const GOLDEN = path.resolve(path.dirname(fileURLToPath(import.meta.url)), 'golden');
const load = f => parseGolden(fs.readFileSync(path.join(GOLDEN, f), 'utf8'));

function expectMatches(actual, file) {
  const mism = findMismatches(actual, load(file));
  expect(mism, `Golden mismatch vs ${file}:\n${mism.join('\n')}`).toEqual([]);
}

test('math golden (solvers + expval) unchanged', () => {
  expectMatches(buildMathGolden(), 'math.golden.json');
});

test('curves golden (figure data for all models/scenarios) unchanged', () => {
  expectMatches(buildCurvesGolden(), 'curves.golden.json');
});

test('fit golden (calculate_lsf recovery) unchanged', () => {
  expectMatches(buildFitGolden(), 'fit.golden.json');
});
