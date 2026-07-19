# bindingsims — Refactor Blueprints

**Status: design only. No implementation.** This document is the master plan the branches will
follow. It is deliberately verbose so that a first-time reader (or a future you) can rebuild the
mental model from scratch. Read [REFACTOR_PLAN.md](REFACTOR_PLAN.md) first for the *why*; this
document is the *how*.

Two rules govern every step:

1. **Behaviour-preserving.** The four existing simulations must produce identical curves,
   populations, pie values, and fit results at every step. The gate is a golden-output test suite
   built in step 0, before any source is touched.
2. **Front-end stays simple and static.** No backend, one HTML page per simulation, deep-linkable
   and embeddable. Machinery (tests, optional build) is dev-time only and never required to *view*
   a sim.

---

## 0. TL;DR of the sequence

Each numbered item is one branch off `main`, merged before the next begins.

| Branch | Goal | Behaviour gate |
|---|---|---|
| `refactor/00-harness-goldens` | test runner + capture golden snapshots of current app | goldens defined |
| `refactor/01-core-math` | move math/expval/solvers into modules, bodies unchanged | goldens pass |
| `refactor/02-render-split` | split `update()` → build / render / table layers | goldens pass |
| `refactor/03-model-registry` | 4 existing models become descriptors; delete `switch(appmode)` | goldens pass |
| `refactor/04-generic-fitter` | pure `fit()` over any model; calcmodes 1–4 identical | fit goldens pass |
| `refactor/05-palette-darkmode` | central palette, CSS vars, theme toggle, theme-aware export | goldens pass (light) |
| `refactor/06-embedding` | drop-in bundle + web component + examples | existing pages unchanged |
| `feature/07-inhibition` | enzyme inhibition (comp/uncomp/noncomp/mixed) | new-model tests |
| `feature/08-unfolding` | 2-state & 3-state unfolding | new-model tests |
| `feature/09-kinetics` | rate curves (time domain) | new-model tests |
| `feature/10-plate-matrix` | 96-well competing-ligand matrix + global fit | new-model tests |

Refactor branches 00–06 change *no observable behaviour*. Feature branches 07–10 only *add*.
The feature order is deliberately low→high effort so the framework is battle-tested on the easy
models before the plate matrix stresses it.

---

## 1. Decision points

All blocking decisions are now **made** (see §1.1). Plain-English explanations of the concepts
behind D1 and D3 are in **§10 (Concepts for non-web-developers)** — read that first if any term
below is unfamiliar. §1.2 lists the smaller best-practice calls.

### 1.1 Decisions made

**D1 — Tooling posture → LEVEL 1 (dev-only node).** Author the code as native **ES modules**
(`import`/`export`). `npm` is used **only on the developer's machine** to run tests (Vitest + jsdom)
and, later, to produce the optional embed bundle. **What ships stays plain static files** — `.htm`,
`.js`, `.css`, no build, no compilation; the files we edit are the files the browser loads. The only
cost: native ES modules must be *served over HTTP* rather than opened by double-clicking the file
(see §10.2), so local preview uses a one-line static server (`python3 -m http.server`); every real
deployment — GitHub Pages, a lab server, an embed — is already HTTP and needs nothing extra.
(Levels 0 "zero tooling" and 2 "full Vite/TypeScript build" were the alternatives; 0 makes the
golden test-gate manual and limits embedding to iframes, 2 adds a build step between edit and view.)

**D2 — Typing → JSDoc-typed JS.** Plain JavaScript annotated with JSDoc comments and type-checked by
`tsc --checkJs` with **no compile step**. Catches the whole "magic index / wrong species order"
class of bugs for free. (Full TypeScript was the alternative; only worth it at D1 Level 2.)

**D3 — Embedding → iframe now, web component at step 06.** The current static site can already be
embedded in any other page via an `<iframe>` today (zero code). The nicer, native-feeling
`<binding-sim model="ligand">` **web component** is built in `refactor/06-embedding`, once the module
refactor has removed the global-name collisions that make drop-in embedding unsafe. Both are
documented; iframe is the interim, the component is the destination. See §7 for the mechanics and
§10.3 for what "iframe" and "web component" mean.

**D4 — New-model data input → deferred to step 10 (decided in principle).** The current fitter reads
a 2-column `x y` textarea, which suits the 1-D models (incl. kinetics). The **96-well matrix** will
add a paste-a-matrix / CSV input plus a signal-model (scale/offset). No action until
`feature/10-plate-matrix`; flagged so it is not a surprise.

### 1.2 Decided — FYI (say the word to change any)

- **Multi-page static site stays.** One `.htm` per sim (great for deep links + embedding). Not a
  single-page app.
- **Keep `?ext` power-user mode** and the `newdecadescale` flag (default `true`). Both paths are
  covered by tests rather than deleted, so nothing silently changes.
- **Golden strictness:** *exact* equality where code is moved verbatim; relative tolerance `1e-12`
  only where float expressions are legitimately re-associated. Any diff beyond that fails the gate.
- **Slider contract frozen.** `expval` + `decadetable` + `expparams` move byte-identical and become
  the single source of truth for slider→value; every conversion routes through them.
- **The uniform solver contract stays** `{ d:[…], m:[…] }`, extended with an optional `signal`
  field for signal-based models (unfolding/kinetics/plate). No breaking change to existing solvers.
- **Feature order:** inhibition → unfolding → kinetics → plate (your call, adopted).
- **Colours become palette tokens** (not RGB literals) in step 05; dark mode rides on that.
- **Each feature branch ships its model + tests + one worked example + a doc note** — no "code now,
  document later".

---

## 2. Target architecture (the shape we are building toward)

Author-time module layout (Level 1). At Level 0 these are the same files as namespaced globals
(`BS.core.expval` etc.) instead of `import`.

```
core/
  expval.js     expval(), decadetable, expparams          — slider↔value contract (frozen numerics)
  solvers.js    solve_cubic_newton/_bisection, root/ODE helpers for future models
  fit.js        fit(model, state, data, opts) → FitResult  — pure, DOM-free
  figure.js     buildFigure(model, state) → FigureSpec     — pure data, no DOM
  render.js     renderFigure(spec, palette, svgRefs); syncDataTable(...); drawPlate(...)
  palette.js    theme tokens → resolved colours (enables dark mode + themed export)
  format.js     number/exponential/superscript formatting, render_text_svg (moved verbatim)
models/
  _template.js          annotated skeleton for new models
  ligand.js  homodimer.js  competing-ligands.js  competing-receptors.js
  (later) inhibition.js  unfolding.js  kinetics.js  competition-matrix.js
  index.js              the registry: { ligand, homodimer, … }
app.js           per-page glue: pick model by id → readState → buildFigure → renderFigure
pages/*.htm      one thin page per sim (loads app.js with a model id)
```

**Data flow (identical for every model, present and future):**

```
DOM sliders ──readState──▶ state ──buildFigure──▶ FigureSpec ──renderFigure──▶ SVG
                              │                        │
                              └────── fit() ◀── data ──┘ (writes recovered params back to sliders)
```

The three `switch(appmode)` blocks in today's code
([update.js curves](../update.js), [calculations.js fun/datalabels](../calculations.js)) collapse
into: *look up the model in the registry, call its `solve`*.

---

## 3. The model descriptor — full contract

This is the heart of the design: **everything the current switches hard-code becomes declarative
data on a descriptor.** All four future models were used to pressure-test the fields below.

```js
/**
 * @typedef {Object} Param
 * @property {string}  key           Canonical name: "K_D", "E_0", "kcat", "dG1"…
 * @property {string}  label         Rich-text token string (see format.js) e.g. "KD"
 * @property {string}  unit          "mol l⁻¹", "s⁻¹", "kJ mol⁻¹", ""
 * @property {[number,number,number]} exp   expparams row [base, minExp, maxExp] for expval()
 * @property {number}  sliderDefault Default raw slider position 0..240
 * @property {boolean} fittable      May the fitter vary it?
 * @property {"continuous"|"discrete"} [kind]   default "continuous"
 * @property {string[]} [choices]    for discrete params (e.g. inhibition mechanism)
 * @property {boolean} [derivedOnly] display-only read-out (ΔG, K_A); no slider
 * @property {(p:Object)=>number} [derive]      compute a derivedOnly value from other params
 */

/**
 * @typedef {Object} Species
 * @property {string}  key           "PL", "P", "f_U"…
 * @property {string}  label         rich-text token string
 * @property {string}  colour        PALETTE TOKEN ("series-1"), never an RGB literal
 * @property {number}  m             stoichiometric multiplicity (as today's m[])
 * @property {boolean} [inPie=true]
 * @property {boolean} [inTable=true]
 */

/**
 * @typedef {Object} XMode          one entry per horizontal-axis choice (today's radio3/4/11)
 * @property {string}  id            "total" | "free" | "specificity" | "time" | "denaturant"
 * @property {string}  label         axis title text
 * @property {"lin"|"log"|"time"} axis
 * @property {string}  independent   which quantity the sweep varies
 */

/**
 * @typedef {Object} YMode          vertical-axis choice (today's radio1/2 and 5/6/7)
 * @property {"absolute"|"relative"|"specificity"|"signal"} id
 * @property {"lin"|"log"} axis
 */

/**
 * @typedef {Object} Marker         dashed vertical reference lines (today's curves[2..8])
 * @property {string}  fromParam     draw a vertical line at this param's value ("K_D", "E_0"…)
 * @property {string}  label
 * @property {string}  colour        palette token
 */

/**
 * @typedef {Object} SolveResult
 * @property {number[]} d            species concentrations, ORDER MATCHES species[]
 * @property {number[]} m            multiplicities (mirror of species[].m; kept for parity)
 * @property {number}  [signal]      optional observable (baseline-modified) for signal models
 */

/**
 * @typedef {Object} Model
 * @property {string}   id
 * @property {string}   title
 * @property {string}   scheme        "P + L ⇌ PL"
 * @property {Param[]}   params
 * @property {Species[]} species
 * @property {XMode[]}   xModes
 * @property {YMode[]}   yModes
 * @property {Marker[]}  [markers]
 * @property {"curve"|"family"|"plate"} layout   default "curve"
 * @property {FamilySpec} [family]     for curve families (one curve per value of a param)
 * @property {PlateSpec}  [plate]      for 2-D plate layout (rows/cols)
 * @property {(p:Object, x:number, ctx:SolveCtx)=>SolveResult} solve   ← the ONLY math
 */
```

### 3.1 Worked example: the existing ligand model (bodies verbatim)

```js
export const ligand = {
  id: "ligand", title: "Ligand binding simulation", scheme: "P + L ⇌ PL",
  params: [
    { key:"E_0", label:"c_P", unit:"mol l⁻¹", exp:[10,-9,-1], sliderDefault:150, fittable:true },
    { key:"S_0", label:"c_L", unit:"mol l⁻¹", exp:[10,-9,-1], sliderDefault:120, fittable:true },
    { key:"K_D", label:"K_D", unit:"mol l⁻¹", exp:[10,-9,-1], sliderDefault: 90, fittable:true },
    { key:"dG",  label:"ΔG",  unit:"kJ mol⁻¹", derivedOnly:true, derive:p=>8.31446*298.15*Math.log(p.K_D)*1e-3 },
    { key:"K_A", label:"K_A", unit:"l mol⁻¹",  derivedOnly:true, derive:p=>1/p.K_D },
  ],
  species: [
    { key:"PL", label:"[PL]", colour:"series-1", m:1 },
    { key:"P",  label:"[P]",  colour:"series-2", m:1 },
  ],
  xModes: [
    { id:"total", label:"Total ligand concentration", axis:"log", independent:"S_0" },
    { id:"free",  label:"Free ligand concentration",  axis:"lin", independent:"S_free" },
  ],
  yModes: [ {id:"absolute",axis:"log"}, {id:"relative",axis:"lin"} ],
  markers: [ {fromParam:"S_0",label:"c_L",colour:"series-5"},
             {fromParam:"E_0",label:"c_P",colour:"series-6"},
             {fromParam:"K_D",label:"K_D",colour:"series-6"} ],
  solve(p, x, ctx) {
    return ctx.xMode === "free"
      ? calculate_ligand_free (p.E_0, x, p.K_D)   // ← today's function, moved unchanged
      : calculate_ligand_total(p.E_0, x, p.K_D);
  },
};
```

The `solve` bodies are literally today's `calculate_*` functions relocated into `core/solvers` and
called from the descriptor. That is what keeps the four models numerically identical — **only the
dispatch changes, never the arithmetic.**

### 3.2 Forward-designed seams (built now, exercised by future models)

These generalizations are cheap to design during the refactor and painful to retrofit, so they go
into the contract in steps 02–03 even though nothing uses them until the feature branches:

- **`family` layout** — draw one curve per value of a chosen param. Needed by *inhibition*
  (a curve per `[I]`) and reused by the *plate* column view. `buildFigure` loops the family param.
- **`signal` yMode + optional `SolveResult.signal`** — an observable with sloping baselines rather
  than a raw concentration. Needed by *unfolding* (CD/fluor. signal) and *kinetics*, *plate*.
- **`time` x-axis** — wire the already-defined but unused `axistype_time`
  ([simulation.js:58](../simulation.js)) so *kinetics* needs no axis surgery.
- **`plate` layout + `drawPlate`** — a render seam (a no-op stub in step 02) that the *matrix*
  model fills with an 8×12 heatmap. The rest of `renderFigure` never learns about plates.
- **Vector-based fit** — `fit()` takes `data` and a `predict(params) → number[]` and minimizes over
  an opaque residual vector, so the *plate* just supplies a length-96 vector; the fitter is unchanged.

---

## 4. Layer contracts (what each core module promises)

**`core/expval.js`** — pure. `expval(raw, base, minExp, maxExp)` and `expval_wrap(raw, exp[])`
return the exact value today's code returns; `decadetable`/`expparams` exported read-only. Frozen by
golden tests; never edited for behaviour.

**`core/solvers.js`** — pure. Houses the moved binding math and the cubic solvers, plus (later)
generic root/exponential helpers. No DOM, no globals. Every solver keeps the `{d,m}` shape.

**`core/figure.js` — `buildFigure(model, state) → FigureSpec`** — pure, no DOM. `state` is
`{ params:{K_D,E_0,…}, xMode, yMode, decadeshift, masses }`. Returns everything render needs:

```js
FigureSpec = {
  curves:   [{ points:[{x,y}], colour, label, dashed, legendSlot|null }],  // main + markers
  pie:      [{ value, colour }],
  populations:[{ key,label,conc,gramConc,percent,colour }],   // the data table rows
  axes:     { xmin,xmax,ymin,ymax, xType, yType, xLabel, yLabel, xTicks, yTicks },
  plate?:   { rows, cols, cells:[[signal]], colourScale },    // only when layout==="plate"
}
```

This is the contract that guarantees rendering is model-agnostic: it must be able to express
*everything* today's `update()` draws (main curves, dashed Kd/conc marker lines, datapoint crosses,
pie, legend-rect vs inline-label, lin/log ticks, white masks, axis titles).

**`core/render.js` — `renderFigure(spec, palette, svgRefs)`** — the only module that touches the SVG
DOM. Knows nothing about binding. Reuses today's SVG element-reuse pattern (create-once, update
attributes) so exported markup stays stable. `palette` supplies resolved colours from tokens.

**`core/fit.js` — `fit(model, state, data, opts) → { values, covariance, residuals, status }`** —
pure. Reuses the existing grid-search (calcmodes 1–3) and Jacobian/Gauss-Jordan covariance
(calcmode 4) *verbatim*, but parameterized by the model's `fittable` params + `solve` instead of the
DOM and `fun_*` closures. The DOM read/write (slider positions in, results out) lives in `app.js`,
not in `fit()`.

**`app.js`** — the only file allowed to read/write the DOM for state. `readState` pulls slider raw
values → `expval` → params; wires radios to xMode/yMode; `writeState` pushes fitted params back to
slider positions (using the inverse of `expval`, as today). Keeps `update()`'s orchestration role
but delegates all computation and rendering.

---

## 5. Step-by-step blueprint (per branch)

Every branch: (1) do the work, (2) run goldens, (3) update the relevant doc, (4) PR to `main`.

### `refactor/00-harness-goldens`
- Add `package.json` + chosen runner (Level 1: Vitest + jsdom). No change to any sim file.
- Load the **current, unmodified** `simulation.js`/`calculations.js`/`update.js` into jsdom against a
  fixture page; drive each sim across the full matrix of {model} × {yMode} × {xMode} ×
  {a few decadeshifts} × {representative slider sets}.
- Snapshot: `curves[]` arrays, data-table populations, `pie[]`, and `calculate_lsf` results for a
  fixed synthetic dataset per model. These files are the **golden reference** for all later steps.
- Also snapshot `expval` over a sweep of raw positions (both `newdecadescale` values).
- Deliverable: `test/golden/*.json` + a `npm test` that re-derives and compares. **No behaviour change.**

### `refactor/01-core-math`
- Create `core/expval.js`, `core/solvers.js`, `core/format.js`; move the bodies **unchanged**.
- At Level 0: wrap in a `BS` namespace; at Level 1: `export`/`import`.
- Pages still load and work. Run goldens (must be exact).

### `refactor/02-render-split`
- Extract `buildFigure` and `renderFigure` + `syncDataTable` + `drawPlate` (stub) from `update()`.
- `update()` becomes: `readState → buildFigure → renderFigure`. Introduce the `FigureSpec` type.
- Wire `axistype_time` and the `signal` yMode plumbing (unused but present). Add `family` loop
  (degenerate: one curve) — no visible effect yet.
- Run goldens (exact; this is the riskiest step — the render seam must reproduce byte-identical SVG
  structure, verified by the DOM smoke snapshot).

### `refactor/03-model-registry`
- Author `models/{ligand,homodimer,competing-ligands,competing-receptors}.js` + `models/index.js`.
- Replace the three `switch(appmode)` blocks with registry lookups. `appmode` integers may remain as
  legacy ids behind the model `id` for one release, then retire.
- Add the **model-contract test** (§6) that runs generic invariants over every registered model.
- Run goldens. **Checkpoint: the behaviour-preserving refactor is complete.**

### `refactor/04-generic-fitter`
- Move `calculate_lsf` internals into `core/fit.js`, parameterized by model. DOM read/write moves to
  `app.js`. Keep calcmodes 1–4 arithmetic identical.
- Run fit goldens for all four models and all calcmodes (exact / 1e-12).

### `refactor/05-palette-darkmode`
- Introduce `core/palette.js` + CSS custom properties; replace `colour1..7` literals and the
  `fill:white` masks/background with tokens. Add a theme toggle (`data-theme`, `localStorage`).
- Make SVG export inline **resolved** colours (clone + computed-style + serialize) so downloads are
  self-contained in either theme; keep a light-theme export default for publication.
- Run goldens in light theme (must match); add dark-theme render smoke tests.

### `refactor/06-embedding`
- Package a self-contained bundle per sim and a `<binding-sim model="…">` web component (Shadow DOM
  for style isolation). Add `examples/iframe.html` and `examples/web-component.html`.
- Write `docs/EMBEDDING.md`. Existing pages unchanged (goldens still pass).

### Feature branches (framework now stable)
Each adds one `models/*.js`, one `pages/*.htm` (from the page template), tests, and a doc + example.

- **`feature/07-inhibition`** — finish the stranded stub. Declare `appmode`-free model with params
  `kcat,Km,E_0,I_0,Ki,Ki'` + discrete `mechanism`. Implement the corrected α/α′ forms (competitive,
  uncompetitive, non-competitive = mixed with Ki=Ki′, mixed) — the current `alpha()`
  ([calculations.js:821](../calculations.js)) conflates these. `layout:"family"` → a curve per
  `[I]`. Cheng–Prusoff IC50↔Ki helper surfaced.
- **`feature/08-unfolding`** — `states:2|3`. `solve(p,D)` returns populations `f_N,(f_I),f_U` in
  `d[]` and a baseline-modified `signal`. `yMode:"signal"`, `xMode:"denaturant"` (LEM) with a thermal
  variant later. Populations flow through the existing table/pie path. Fit ΔG(s), m-value(s),
  baselines; report Cm/Tm.
- **`feature/09-kinetics`** — `xMode:"time"` (seam already wired). `solve(p,t)` = approach-to-eq
  exponential + baseline; optional assoc/dissoc phases; k_obs-vs-[L] replot helper to extract
  k_on/k_off. Fit rate constants + amplitudes.
- **`feature/10-plate-matrix`** — `layout:"plate"`, `plate:{rows:[L…],cols:[L'…]}`. Per-well `solve`
  reuses `calculate_ligands` (the cubic already exists). `drawPlate` renders the 8×12 heatmap;
  column view uses `family`. Global fit via the vector-based `fit()` (residual vector length = 96).
  New matrix/CSV data input (D4).

---

## 6. Testing strategy

Runner per D1 (Level 1: Vitest + jsdom; Level 0: in-browser QUnit page). Layers:

- **Golden / regression (step 00, the reproducibility gate).** Snapshots of curves, populations,
  pie, and fit results per model × mode combination; must not drift through steps 01–06.
- **Model-contract suite (parameterized over the registry).** For *every* registered model, assert:
  species order length == `solve().d.length`; unique param keys; valid `exp` rows; mass/really
  population conservation on random inputs; markers reference real params; `fittable` params round-trip
  through `expval`. New models inherit this coverage automatically — the main guardrail against
  future bloat/regressions.
- **Unit tests** (from REFACTOR_PLAN §7): solvers (conservation, limiting cases, cubic
  Newton/bisection agreement, stable vs unstable quadratic form), `expval` exactness + monotonicity +
  decadeshift, `fit` parameter recovery + covariance sanity + guard rails, `sum_squared_residuals`,
  `tinv_95`.
- **Render smoke (jsdom).** Each page's `renderFigure` runs without throwing and emits the expected
  node counts — cheap protection for the step-02 render split and step-05 palette work.
- **Per-new-model tests** (REFACTOR_PLAN §7E): inhibition reduces to MM at I=0 and Cheng–Prusoff
  round-trips; unfolding f_U=½ at Cm and 3-state→2-state limit; kinetics hits A_eq/t½ and k_obs
  replot slopes; plate per-well matches `calculate_ligands` and global fit recovers K_D,K_D′.

**Bloat control:** the contract suite + goldens mean a new model can't merge unless it satisfies the
invariants and doesn't perturb existing snapshots. Keep tests colocated (`test/` mirroring `core/`
and `models/`), one golden file per model.

---

## 7. Static site & embedding (feature request #2)

**Recommendation: stay static — it is the right architecture here.** These are pure client-side
simulations: no server state, no auth, no database. Static means free hosting (GitHub Pages),
trivial caching, deep-linkable pages, and — once the module refactor removes the global-namespace
collisions — safe embedding. The only things static can't do (server-side compute, DB, SEO of
dynamic content) are irrelevant to a binding simulator. I would not move off static.

Embedding options, cheapest → cleanest:

1. **iframe (available now).** Host the site; others embed `<iframe src=".../ligand.htm">`. Zero
   code, fully isolated CSS/JS. Downside: manual sizing, feels bolted-on. **Documented as the
   baseline in step 06.**
2. **Folder drop-in.** Copy a sim's files into a subfolder of another static site. Works because
   paths are relative — but only *safe* after the refactor namespaces/module-scopes the globals
   (today's `S_0`, `update`, etc. would collide with a host page). Enabled by steps 01–03.
3. **Web component `<binding-sim model="ligand" params="K_D=1e-6">` (step 06).** Custom element with
   Shadow DOM for style isolation; one `<script>` include, then drop the tag anywhere. The cleanest
   "embed in other websites" story and the reason step 06 produces a bundle. Needs D1 ≥ Level 1.
4. **npm package** — deferred; only if third-party app builders want it.

So: static stays; iframe now; web-component as the aspirational embed built on the modularized core.
This is why the refactor (killing global collisions) is a prerequisite for *good* embedding, not
just cleanliness.

---

## 8. Documentation & examples for future extenders

Shipped alongside the code so the repo teaches its own extension:

- **`docs/ARCHITECTURE.md`** — the module map + the layer contracts from §4 (distilled from this
  doc), kept current as the refactor lands.
- **`docs/ADDING_A_MODEL.md`** — a step-by-step tutorial: copy `models/_template.js`, fill the
  descriptor, add a page from `pages/_template.htm`, add a golden test, register it. Written against
  the inhibition model as the concrete example.
- **`models/_template.js`** — a heavily-commented skeleton descriptor (every field explained).
- **Self-documenting model files** — each `models/*.js` opens with the biophysics, the equations,
  and a literature reference, so the science is legible next to the code.
- **`examples/`** — `iframe.html` and `web-component.html` working embeds.
- **`docs/EMBEDDING.md`** — how to host and embed (from §7).
- The existing GPL headers and the ProtSim citation are preserved in every moved file.

---

## 9. Risks & how the plan contains them

- **Silent numeric drift** → the golden gate on every step; exact equality for verbatim moves.
- **Render seam (step 02) diverging from current SVG** → DOM smoke snapshots + goldens on the drawn
  coordinates, not just the model arrays.
- **Fitter behaviour change (step 04)** → fit goldens across all calcmodes before/after.
- **Over-engineering** → the front end stays multi-page static; machinery is dev-only; the contract
  suite blocks bloat; new models are *data* (descriptors), not new control flow.
- **Scope creep in the plate model** → it is last on purpose, after the framework has absorbed
  families, signal baselines, and vector fitting from the three easier models.

---

## 10. Concepts for non-web-developers

Plain-English background for the terms used above, aimed at a scientist who is not a web developer.
Nothing here changes the plan; it exists so this document can remind a future reader what the jargon
means.

### 10.1 How the code loads today, and why globals are a problem

Each page includes scripts with tags like `<script src="calculations.js"></script>`. These are
*classic scripts*: every `var` and `function` written at the top of the file becomes a **global** — a
name attached to the browser's shared `window` object and visible to every other script on the page.
That is *why* the app works today (`update.js` can call `calculate_ligand_total` from
`calculations.js` because both are globals), and also why it is brittle: any two names can collide,
and dropping these scripts into someone else's web page risks clashing with *their* code (their page
might also have a variable called `update` or `S_0`). Removing this global sharing is what makes the
code safe to embed and easier to test.

### 10.2 ES modules, and "served over HTTP"

- **ES module** — "ES" is *ECMAScript*, the official name for JavaScript. An ES module is simply a
  `.js` file that uses two keywords: `export` ("here is what this file shares") and `import` ("here is
  what this file needs"), e.g. `export function expval(){…}` and
  `import { expval } from './core/expval.js'`. Nothing leaks into the global pool unless explicitly
  exported. A page opts in with `<script type="module" src="app.js">`. Same JavaScript language — the
  difference is that dependencies are stated out loud instead of assumed through globals. Benefits:
  no name collisions (safe embedding), a visible dependency graph, and the ability for a test tool to
  import one function and test it in isolation.

- **"Served over HTTP"** — When you double-click an HTML file, the browser opens it with a `file://`
  address. Classic scripts work over `file://`; **ES modules do not** — browsers block modules loaded
  from `file://` for security reasons. "Served over HTTP" just means the files are delivered by a web
  server at an `http://`/`https://` address. Two everyday cases, both trivial and both still *just
  static files*:
  - **Deploy:** GitHub Pages (`https://…github.io/bindingsims/`), a lab server, etc. — already HTTP,
    so modules work with zero extra effort.
  - **Local preview:** run one command in the folder — `python3 -m http.server` — then open
    `http://localhost:8000/ligand.htm`. That "server" only hands out files; it does no computation.
  - Files may live in subfolders (`core/`, `models/`); imports use relative paths. "Same origin"
    means "same website," not "same folder." The only thing lost versus today is *double-click to
    preview*; the one-line local server replaces it.

- **What "Level 1 tooling" ships** — exactly what ships today: plain `.htm`, `.js`, `.css`, no build,
  no compilation. `npm`/node is used **only on the developer's machine** to run the test suite; it
  never touches what a visitor downloads.

### 10.3 Embedding: iframe vs web component

"Embedding" means making a simulation appear *inside another web page* (a lab site, a course page, a
blog post).

- **iframe (available now).** An `<iframe>` is an HTML tag that embeds one web page inside another,
  like picture-in-picture. Host bindingsims somewhere, and the other page adds one line:
  `<iframe src="https://…/ligand.htm" width="700" height="900"></iframe>`. The sim appears in a boxed
  sub-window. This works with the *current* static site, no code changes. Downsides: it is a
  fixed-size rectangle you size by hand, it reads as a bolted-on window rather than part of the page,
  and host↔sim communication is clumsy.

- **Web component (the "clean" embed, step 06).** A web component lets us define our *own* HTML tag.
  After including one script, the host page writes `<binding-sim model="ligand" kd="1e-6"></binding-sim>`
  and the sim renders inline like a native element. "Clean" means three concrete things: (a) it sits
  *in* the page rather than in a boxed iframe; (b) it is configurable through attributes
  (`model=`, `kd=`); (c) it is shielded from the host page's styles by a browser feature called
  **Shadow DOM** — a private styling bubble, so the host's CSS cannot accidentally restyle the sim and
  the sim's CSS cannot leak out. It needs the module refactor first (no global collisions) plus small
  packaging, which is why it lands in `refactor/06-embedding` rather than now.

### 10.4 The test vocabulary used above

- **Golden / snapshot test** — record the current program's exact outputs (the curve arrays,
  populations, fit results) to files once, then on every later change re-run and compare against those
  recorded files. Any difference fails the test. This is how "numerically identical" is enforced.
- **Unit test** — checks one function in isolation (e.g. "the cubic solver returns the right root").
- **Contract test** — one test that runs the same invariants over *every* registered model (species
  count matches solver output, populations conserve mass, etc.), so new models get baseline coverage
  automatically.
- **jsdom** — a fake browser environment that runs in the terminal, so tests can exercise the page
  logic (including SVG DOM) without opening a real browser.
- **Vitest** — the test runner (the program that finds and executes the tests and reports pass/fail).
- **CI** — "continuous integration": running the test suite automatically on every push (e.g. via
  GitHub Actions) so a regression is caught without anyone remembering to run tests by hand.
