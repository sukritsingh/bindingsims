// Loads an unmodified simulation page into jsdom and runs its scripts exactly as
// a browser would, so goldens capture real behaviour rather than a reimplementation.
//
// The three page scripts (simulation.js / calculations.js / update.js) are still
// CLASSIC scripts: injected as inline <script> elements so their top-level
// `var`/`function` declarations attach to `window` as globals — exactly how the
// app works in the browser today.
//
// In refactor/01-core-math the pure math/format/expval helpers moved into the
// `core/` ES modules. In the browser, `core/bootstrap.js` imports those modules
// and re-exposes their functions on `window` (the compatibility shim), while the
// moved functions read shared state (`newdecadescale`, `decadeshift`, `svg_xmlns`,
// `delete_all_children`, `appmode`) off `window`. jsdom does not execute
// `type="module"` scripts with filesystem imports, so THIS loader reproduces the
// shim in Node instead: it imports the same core modules, bridges the shared
// state from Node's global scope to the active jsdom window (so the moved
// functions resolve it live), and copies the module exports onto the jsdom
// window. The committed golden JSON stays frozen — only this loader changed.

import { JSDOM } from 'jsdom';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

import * as expvalMod from '../../core/expval.js';
import * as solversMod from '../../core/solvers.js';
import * as formatMod from '../../core/format.js';
import { models, modelByAppmode, datalabels } from '../../models/index.js';

export const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..', '..');

/** Read a file relative to the repository root. */
export function readRepo(rel) {
  return fs.readFileSync(path.join(ROOT, rel), 'utf8');
}

// Which pages map to which appmode global. index.htm is list-only (no figure).
export const PAGES = {
  ligand:     { file: 'ligand.htm',     appmode: 'appmode_ligand' },
  homodimer:  { file: 'homodimer.htm',  appmode: 'appmode_homodimer' },
  ligands:    { file: 'ligands.htm',    appmode: 'appmode_ligands' },
  receptors:  { file: 'receptors.htm',  appmode: 'appmode_receptors' },
};

const SCRIPTS = ['simulation.js', 'calculations.js', 'update.js'];

// The shared state that the moved core functions read as free globals. In the
// browser these resolve to `window`; here we bridge Node's global scope to the
// active jsdom window with live getters, so mutations made through the window
// (e.g. tests setting `w.decadeshift`) are seen by the imported modules.
const BRIDGED_GLOBALS = [
  'newdecadescale', 'decadeshift',        // expval.js
  'appmode', 'appmode_ligands',           // solvers.js
  'svg_xmlns', 'delete_all_children',     // format.js
  'document',                             // format.js (render_text_svg)
];

let activeWindow = null;

for (const name of BRIDGED_GLOBALS) {
  Object.defineProperty(globalThis, name, {
    configurable: true,
    get() { return activeWindow ? activeWindow[name] : undefined; },
  });
}

/**
 * Load a page and return its live `window` with the sim globals populated and
 * init() already run. `opts.ext` appends ?ext to enable power-user mode.
 */
export function loadSim(pageKey, opts = {}) {
  const page = PAGES[pageKey];
  if (!page) throw new Error(`unknown page: ${pageKey}`);

  let html = readRepo(page.file);
  // Remove the external <script src> tags (classic and the module bootstrap) and
  // the inline preinit call so we control load order and initialisation ourselves.
  html = html.replace(/<script src="[^"]*"><\/script>\s*/g, '');
  html = html.replace(/<script type="module" src="[^"]*"><\/script>\s*/g, '');
  html = html.replace(/<script>\s*preinit\([^)]*\);\s*<\/script>/g, '');

  const url = 'http://localhost/' + page.file + (opts.ext ? '?ext' : '');
  const dom = new JSDOM(html, { runScripts: 'dangerously', pretendToBeVisual: true, url });
  const { window } = dom;
  window.alert = () => {};          // calcmode 4 alerts its report; swallow it

  // Inject the classic page scripts; their globals populate this window.
  for (const f of SCRIPTS) {
    const s = window.document.createElement('script');
    s.textContent = readRepo(f);
    window.document.body.appendChild(s);
  }

  // Reproduce core/bootstrap.js: bridge shared state to this window and re-expose
  // the moved functions on it, so the classic code can call them by bare name.
  activeWindow = window;
  Object.assign(window, expvalMod, solversMod, formatMod);

  // Model registry (refactor/03), mirroring core/bootstrap.js: the classic code
  // reads modelByAppmode(...) and datalabels[] off the window.
  window.models = models;
  window.modelByAppmode = modelByAppmode;
  window.datalabels = datalabels;

  window.preinit(window[page.appmode]);
  window.init();
  return window;
}
