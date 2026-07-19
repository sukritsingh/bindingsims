// Loads an unmodified simulation page into jsdom and runs its scripts exactly as
// a browser would, so goldens capture real behaviour rather than a reimplementation.
//
// The three sim scripts are injected as inline <script> elements (classic scripts),
// so their top-level `var`/`function` declarations attach to `window` as globals —
// which is precisely how the app works in the browser today. `preinit()`/`init()`
// are then called by hand (we strip the page's inline preinit + rely on no onload).
//
// When the refactor moves code into modules, only THIS loader changes to match the
// new structure; the committed golden JSON stays frozen. That is the whole point.

import { JSDOM } from 'jsdom';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

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

/**
 * Load a page and return its live `window` with the sim globals populated and
 * init() already run. `opts.ext` appends ?ext to enable power-user mode.
 */
export function loadSim(pageKey, opts = {}) {
  const page = PAGES[pageKey];
  if (!page) throw new Error(`unknown page: ${pageKey}`);

  let html = readRepo(page.file);
  // Remove the external <script src> tags and the inline preinit call so we control
  // load order and initialisation ourselves.
  html = html.replace(/<script src="[^"]*"><\/script>\s*/g, '');
  html = html.replace(/<script>\s*preinit\([^)]*\);\s*<\/script>/g, '');

  const url = 'http://localhost/' + page.file + (opts.ext ? '?ext' : '');
  const dom = new JSDOM(html, { runScripts: 'dangerously', pretendToBeVisual: true, url });
  const { window } = dom;
  window.alert = () => {};          // calcmode 4 alerts its report; swallow it

  for (const f of SCRIPTS) {
    const s = window.document.createElement('script');
    s.textContent = readRepo(f);
    window.document.body.appendChild(s);
  }

  window.preinit(window[page.appmode]);
  window.init();
  return window;
}
