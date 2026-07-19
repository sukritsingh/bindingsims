// Regenerates the committed golden files from the CURRENT code.
//
//   npm run golden
//
// Run this only when you have deliberately, correctly changed behaviour and want
// the goldens to reflect the new truth. During the refactor the goldens must NOT
// be regenerated — the whole point is that the tests keep comparing against the
// behaviour captured here in step 00.

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { buildCurvesGolden, stableStringify, compactStringify } from '../helpers/scenarios.js';
import { buildMathGolden, buildFitGolden } from '../helpers/mathCases.js';

const GOLDEN_DIR = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..', 'golden');

function write(name, text) {
  const file = path.join(GOLDEN_DIR, name);
  fs.writeFileSync(file, text + '\n');
  console.log('wrote', path.relative(process.cwd(), file));
}

write('math.golden.json', stableStringify(buildMathGolden()));
write('curves.golden.json', compactStringify(buildCurvesGolden()));   // compact: large
write('fit.golden.json', stableStringify(buildFitGolden()));
console.log('done.');
