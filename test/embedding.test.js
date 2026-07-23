// Smoke tests for the <binding-sim> web component (refactor/06-embedding).
// The component wraps an isolated iframe pointing at a sim page; we verify it
// registers, resolves each model to its page, honours ?ext, and errors clearly
// on an unknown model. (The sim pages themselves are covered by the goldens.)

import { test, expect } from 'vitest';
import { JSDOM } from 'jsdom';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const src = fs.readFileSync(path.join(ROOT, 'binding-sim.js'), 'utf8');

// Fresh jsdom (own custom-element registry) with binding-sim.js loaded as an
// inline classic script, mirroring how a host page includes it.
function makeWindow() {
  const dom = new JSDOM('<!DOCTYPE html><body></body>', {
    runScripts: 'dangerously',
    url: 'http://localhost/app/',
  });
  const s = dom.window.document.createElement('script');
  s.textContent = src;
  dom.window.document.body.appendChild(s);
  return dom.window;
}

function mount(w, attrs) {
  const el = w.document.createElement('binding-sim');
  for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, v);
  w.document.body.appendChild(el);
  return el;
}

test('registers the <binding-sim> custom element', () => {
  const w = makeWindow();
  expect(typeof w.customElements.get('binding-sim')).toBe('function');
});

test('renders an isolated iframe pointing at the model page', () => {
  const w = makeWindow();
  const el = mount(w, { model: 'ligand' });
  const iframe = el.shadowRoot.querySelector('iframe');
  expect(iframe).toBeTruthy();
  expect(iframe.src).toBe('http://localhost/app/ligand.htm');
  expect(iframe.getAttribute('title')).toContain('ligand');
});

test('all four models resolve to their pages', () => {
  const w = makeWindow();
  for (const m of ['ligand', 'homodimer', 'ligands', 'receptors', 'inhibition']) {
    const el = mount(w, { model: m });
    expect(el.shadowRoot.querySelector('iframe').src).toBe('http://localhost/app/' + m + '.htm');
  }
});

test('ext attribute appends ?ext', () => {
  const w = makeWindow();
  const el = mount(w, { model: 'receptors', ext: '' });
  expect(el.shadowRoot.querySelector('iframe').src).toBe('http://localhost/app/receptors.htm?ext');
});

test('unknown model shows an error and no iframe', () => {
  const w = makeWindow();
  const el = mount(w, { model: 'bogus' });
  expect(el.shadowRoot.querySelector('iframe')).toBeNull();
  expect(el.shadowRoot.querySelector('.error')).toBeTruthy();
});
