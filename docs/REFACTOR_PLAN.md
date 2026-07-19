# bindingsims — refactor & extension plan

Status: planning only. No behaviour on `main` is changed by this document.
Target branch for all work: feature branches off `main` (see §8). This doc lives on
`claude/binding-curve-refactor-5ad60d`.

The guiding constraint throughout: **existing UI, existing models, sliders, Kd fits, and SVG
export must remain numerically identical and behaviourally reproducible.** Everything below is
designed so that the current four simulations produce bit-for-bit the same curves and populations
before and after any refactor. The mechanism for guaranteeing that is a golden-output test suite
captured *before* touching anything (§7).

---

## 1. What the repository is today

A static, build-less web app. Each simulation is one `.htm` page that pulls in three shared
scripts in a fixed order:

```
simulation.js   → globals, slider machinery, DOM/SVG helpers, expval mapping, SVG export
calculations.js → the binding math (pure-ish), cubic solvers, and calculate_lsf() (the fitter)
update.js       → one 1580-line update() that recomputes curves AND renders the whole figure
```

Pages: [index.htm](../index.htm), [ligand.htm](../ligand.htm), [ligands.htm](../ligands.htm),
[receptors.htm](../receptors.htm), [homodimer.htm](../homodimer.htm). Plus two orphans:
[enzyme.htm](../enzyme.htm) (a non-functional stub) and [extensions.htm](../extensions.htm) (docs
for the `?ext` power-user mode).

### The four working models

| appmode | page | math (calculations.js) | form |
|---|---|---|---|
| `appmode_ligand` (2) | ligand.htm | `calculate_ligand_total` / `_free` | quadratic |
| `appmode_homodimer` (3) | homodimer.htm | `calculate_homodimer` / `_free` | quadratic |
| `appmode_ligands` (0) | ligands.htm | `calculate_ligands` (cubic) | 1 protein + 2 competing ligands |
| `appmode_receptors` (1) | receptors.htm | `calculate_ligands` / `_free` | 2 competing proteins + 1 ligand |

Every solver returns the same shape: `{ d: [species concentrations…], m: [stoichiometric
multiplicities…] }`. That uniform contract is the single best design decision in the codebase and
the refactor should keep it verbatim.

### The three core data-flow globals

- **appmode** — an integer selecting the model. It drives a big `switch` in three separate places
  (curve generation in `update.js`, function selection + `datalabels` in `calculations.js`).
- **Slider indices** — parameters are identified by *magic integers* (`3, 5, 7, 9, 10`, with
  `8, 16, 17, 18` for derived read-outs). The same integer keys the HTML element ids
  (`slider5`, `value5`, `value_input5`, `fixval5`, `fixlabel5`), the `expparams` row, the
  `slider_input()` switch, the `valuespan_click()` switch, and the `calculate_lsf()` switch.
- **Physical globals** — `S_0, E_0, K_D, K_D2, Q_0` are mutable module globals in
  [simulation.js:42](../simulation.js) that everything reads and writes.

Mapping of the magic indices (worth pinning down because a refactor must preserve it):

| index | global | meaning | pages |
|---|---|---|---|
| 3 | `S_0` | ligand L′ conc (ligands) / substrate-like | all |
| 5 | `E_0` | total protein P conc | all |
| 7 | `K_D` | first dissociation constant | all |
| 9 | `K_D2` | second dissociation constant | ligands, receptors |
| 10 | `Q_0` | total ligand L conc | ligands, receptors |
| 8 | — | ΔG derived from `K_D` (read-only) | display |
| 16 | — | K_A = 1/`K_D` (read-only) | display |
| 17, 18 | — | ΔG′, K_A′ from `K_D2` (read-only) | display |

---

## 2. Clunkiness & brittleness (what to refactor, and why)

Ordered roughly by how much they block extension.

### 2.1 `update()` is a god-function (highest priority)
[update.js](../update.js) is a single function that, per `appmode`, (a) generates the curve
point arrays, (b) computes the equilibrium data-table numbers, (c) builds the pie chart,
(d) lays out legends, (e) draws axes/ticks/labels, and (f) writes SVG DOM. Model logic and
rendering are fused. **Consequence:** adding a model means adding a ~200-line `case` that repeats
axis/pie/legend boilerplate, and any rendering fix must be duplicated across four cases.

*Improvement:* split into three layers with clean boundaries —
`buildCurves(model, state) → {curves, labels, colours, legends, pie, axes}` (pure data) →
`renderFigure(figureSpec) → SVG` (pure rendering) → `syncDataTable(model, populations)` (DOM).
The model cases shrink to "how do I turn parameters into species", which they already express in
`calculations.js`.

### 2.2 Model identity is scattered across `switch(appmode)` statements
Adding a model today requires touching: `simulation.js` (appmode constant, `expparams`,
`slider_input` cases, `valuespan_click`, `valueinput_*`), `calculations.js` (`datalabels`,
`fun_*` closures, the `switch(appmode)` in `calculate_lsf`), `update.js` (a new `case`), plus a
new HTML page. Five files, many switches, easy to get subtly wrong.

*Improvement:* a **model registry** — one descriptor object per model (see §4). The switches become
lookups keyed by `model.id`.

### 2.3 Magic slider indices
Parameters have no names, only integers, and the integer is duplicated in ≥5 locations per
parameter. Any new parameter needs a free index that does not collide with the derived-value
indices (8, 16, 17, 18). Brittle and non-obvious for an educational codebase people want to read.

*Improvement:* give parameters string keys (`"K_D"`, `"E_0"`) in the descriptor; derive DOM ids
and `expparams` lookups from the key. Keep the integer ids in HTML as a compatibility shim if
needed, but stop hand-maintaining the switches.

### 2.4 The fitter (`calculate_lsf`) is DOM-bound and brute-force
[calculations.js:134](../calculations.js) reads sliders/checkboxes/radios directly from the DOM,
does a nested grid search over integer slider positions (0–240) for calcmodes 1–3, writes results
back into slider positions, and recurses. The genuinely good part — calcmode 4's Jacobian +
Gauss-Jordan covariance/standard-error/t-based CI — is buried, marked "unfinished, may be
unstable", and only reachable in `?ext` mode.

*Improvements (without changing numerics of the default path):*
- Extract a pure `fit(model, params, data, options) → {values, covariance, residuals}` that takes
  a `solve` function and a parameter spec, not the DOM.
- Keep the existing grid-search verbatim as one strategy (so calcmode 1–3 stay identical); expose
  the Jacobian/covariance path as the "report standard errors" step for any model.
- This is what makes the *new* models fittable "for free" — the fitter stops knowing about
  specific `fun_ligand_*` closures.

### 2.5 Hard-coded colours and white masks (this is the dark-mode blocker — see §3)
Curve colours are RGB string literals in `update.js` (`colour1..colour7`,
[update.js:84](../update.js)). The plot background is a `fill:white` `<rect>` in every HTML page.
Three opaque `fill:white` "mask" rectangles ([update.js:1287+](../update.js)) clip curves that
run outside the plot box. Axis lines are `stroke:black`, data-point crosses `fill/stroke:black`,
legend text `fill:white`. None of it is themeable.

### 2.6 `expval` / `decadetable` slider mapping is implicit and untested
[simulation.js:497](../simulation.js) maps slider integer → physical value via a 30-entry decade
table (`newdecadescale=true`) or a power law (`false`). This is the definition of what a slider
position *means* and therefore of reproducibility, yet it is a bare global with no tests and a
dev-only alternate branch.

*Improvement:* leave the function untouched (it is the numeric contract) but wrap it, document it,
and pin it with golden tests (§7). Any refactor routes all slider→value conversions through it.

### 2.7 Dead / half-finished code
- `enzyme.htm` references `appmode_enzyme` (never defined) and undeclared globals
  (`inhib_mech, kcat, E0, Km, I0, Ki`). The enzyme helpers at
  [calculations.js:821](../calculations.js) (`alpha`, `v0`, `mm_rate_curve`, `ic50_to_Ki`) read
  those same undeclared globals. This is a stranded first attempt at the inhibition model (§5.3).
- `extpiemode` / thesis-only pie labels ([update.js:1058](../update.js)) — keep but isolate.

*Improvement:* fold the enzyme stub into the model registry properly (§5.3) or delete it; don't
leave it half-wired.

### 2.8 Smaller sharp edges
- `update()` re-reads `mass1/2/3` and every radio from the DOM on every call; state lives in the
  DOM rather than in a model object.
- Numerical guards are ad hoc: `calculate_ligand_total` documents a stable vs unstable form inline;
  `solve_cubic_newton` bails after 20 cycles and falls back to bisection. Fine, but untested.
- SVG export ([simulation.js:820](../simulation.js)) does string `.replace()` surgery on
  `innerHTML` and opens a blob in a new tab; fragile to any markup change and not theme-aware.
- No package.json, no test runner, no CI.

---

## 3. Dark mode — feasibility and effort

**Verdict: a real but bounded lift. Small–medium (~half a day) *if* colour centralisation
(§2.5) is done first; otherwise it stays whack-a-mole.** The reason it isn't trivial is the
plotting-colour coupling you already anticipated: the figure hard-codes white/black in ~10 places,
and three *opaque white mask rectangles* are load-bearing — they hide curve segments that leave the
plot area, so they must always equal the figure background colour, not just "white".

What a clean implementation looks like:

1. **Central palette.** Replace the `colour1..7` literals and the white/black constants with a
   single palette object read from CSS custom properties:
   `--fig-bg, --fig-axis, --fig-text, --series-1..7, --datapoint`. Light and dark values live in
   `:root` and `:root[data-theme="dark"]` (plus `@media (prefers-color-scheme: dark)`).
2. **Figure background + masks** switch from `fill:white` to `fill: var(--fig-bg)` (mask rects and
   the HTML background `<rect>`). This is the key fix: masks track the theme automatically.
3. **Series colours** must stay distinguishable on dark. The current reds/blues/greens read poorly
   on black — pick a dark-mode variant palette (lighter, desaturated) rather than reusing the same
   seven RGBs. Populations/pie and legend swatches inherit from the same palette, so the data table
   and pie stay consistent for free.
4. **Export.** SVG export must bake the *resolved* colours (getComputedStyle) into the file so a
   downloaded dark-mode figure is self-contained. This is the fiddliest 20%: the current export is
   string surgery on `innerHTML`; it should instead clone the SVG node, inline computed styles, and
   serialise. (A light-background export option is worth keeping as the default for publication.)
5. **Toggle.** A checkbox that stamps `data-theme` on `<html>`; persist in `localStorage`.

Effort ranking: palette + CSS vars (small) → mask/bg swap (small) → dark series palette (design
choice, small) → theme-aware export (medium, the only real work). **Recommendation:** do the
colour centralisation as part of the render-layer split (§2.1/§2.5) and add the toggle afterwards;
don't attempt dark mode against the current inlined colours.

---

## 4. The improved framework (target architecture)

A **model registry** of descriptor objects, one per simulation, consumed by a generic engine.
No framework, no build step required — just split the existing three files into small modules and
add a `models/` folder. (If you later want a build/test toolchain, ES modules + Vitest is the
lightest option; not required for the design.)

```
core/
  expval.js        // expval + decadetable + expparams  (UNCHANGED numerics; the slider contract)
  solvers.js       // solve_cubic_newton/_bisection, future ODE/root helpers
  fit.js           // pure fit(model, paramSpec, data, opts) → {values, covariance, residuals}
  figure.js        // buildCurves(model, state) → figureSpec  (pure data)
  render.js        // renderFigure(figureSpec, palette) → SVG ; syncDataTable(...)
  palette.js       // theme palette (enables §3)
models/
  ligand.js  homodimer.js  competing-ligands.js  competing-receptors.js
  competition-matrix.js  kinetics.js  inhibition.js  unfolding.js
app.js             // wires a page to its model: reads sliders → state → figure → render
```

A **model descriptor** declares everything the switches currently hard-code:

```js
export const ligand = {
  id: "ligand",
  title: "Ligand binding simulation",
  scheme: "P + L ⇌ PL",
  params: [                         // replaces magic indices + expparams rows
    { key: "E_0",  label: "c_P", unit: "mol l⁻¹", exp: [10,-9,-1], slider: 5,  fit: true },
    { key: "S_0",  label: "c_L", unit: "mol l⁻¹", exp: [10,-9,-1], slider: 3,  fit: true },
    { key: "K_D",  label: "K_D", unit: "mol l⁻¹", exp: [10,-9,-1], slider: 7,  fit: true },
  ],
  species: [                        // replaces datalabels + labels + colours + m[]
    { key: "PL", label: "[PL]", colour: "series-1", m: 1 },
    { key: "P",  label: "[P]",  colour: "series-2", m: 1 },
  ],
  xModes: [                         // replaces xscale_alternative branches
    { id: "total", label: "Total ligand concentration", axis: "log"  },
    { id: "free",  label: "Free ligand concentration",  axis: "lin"  },
  ],
  solve(p, x, xMode) {              // the EXISTING calculate_ligand_* body, verbatim
    return xMode === "free"
      ? calculate_ligand_free(p.E_0, x, p.K_D)
      : calculate_ligand_total(p.E_0, x, p.K_D);
  },
};
```

The engine then:
- builds sliders/labels/data-table rows from `params`/`species` (no per-page magic ids to maintain);
- generates curves by sweeping `x` and calling `model.solve` (one loop, not four `case`s);
- fits by handing `model.solve` + the `fit:true` params to `core/fit.js` (no `fun_*` closures);
- renders via `core/render.js` using `palette.js` (themeable).

Crucially, `solve` bodies are the **current** `calculate_*` functions moved unchanged, so the four
existing models stay numerically identical; the registry only replaces the *dispatch*, not the math.

---

## 5. New models — specs, and how to add each in *both* frameworks

For each: the biophysics, the parameters, the fittable quantities, the plot, and the "current
framework" vs "improved framework" cost.

### 5.1 Competing ligands as a 2D matrix (96-well plate)
**Biophysics.** One protein P, two competing ligands L and L′, but instead of a single titration,
a plate: rows = a titration series of L, columns = a titration series of L′ (e.g. 8 rows × 12
columns). Each well is an independent equilibrium. **The core math already exists:**
`calculate_ligands(c_P, c_L, c_L′, K_D, K_D′)` ([calculations.js:87](../calculations.js)) solves
exactly this competition (the cubic). The new work is *sweeping it over a grid* and *fitting the
whole plate at once*.

- **Parameters:** `E_0` (P), `K_D`, `K_D′`, per-row `[L]` values, per-column `[L′]` values, plus a
  signal model (readout = scale·[PL] + offset, or fraction bound). Fit `K_D, K_D′` (and
  signal scale/offset) globally against all 96 readouts simultaneously — this is a global fit, the
  fitter operating on a residual vector of length 96 instead of N.
- **Plot:** two complementary views — (a) the plate as an 8×12 heatmap of predicted signal
  (a new render primitive: a grid of coloured cells with a colour scale), and (b) family-of-curves
  overlays (signal vs [L] for each [L′] column). "How the signal changes across the plate" = the
  heatmap; the pie/data-table generalises to a hovered/selected well.
- **Current framework cost: high.** No 2D concept exists anywhere — `update()` assumes a single
  1D curve sweep, `datapoints` is a flat x,y list, the fitter iterates `datapoints[i].x`. You'd be
  bolting a plate renderer and a matrix data parser onto the god-function. Realistically a new
  appmode plus substantial special-casing in `update()` and `calculate_lsf`.
- **Improved framework cost: medium.** Add `models/competition-matrix.js` whose `solve` returns the
  per-well signal from `calculate_ligands`, declare a `layout: "plate"` so `render.js` draws the
  heatmap, and let `core/fit.js` (already residual-vector based) do the global fit. Data input
  becomes a pasted matrix; the generic fitter needs no change because it already sums squared
  residuals over an arbitrary vector.
- **This is the single strongest argument for doing the refactor first** — it is painful in the
  current design and natural in the target one.

### 5.2 Kinetic curves (fitting rates)
**Biophysics.** Time-domain rather than equilibrium. Cover at least: (a) pseudo-first-order
approach to equilibrium, signal(t) = A_eq·(1 − e^(−k_obs·t)) + baseline; (b) association +
dissociation phases (k_obs = k_on·[L] + k_off, plus a dissociation-only decay); optionally (c)
enzyme progress curves. Fit `k_on, k_off` (or `k_obs`, `A_eq`) and baselines.

- **Parameters:** rate constants, amplitudes, baselines; x-axis is time.
- **Infrastructure already half-present:** `axistype_time = 2` is defined
  ([simulation.js:58](../simulation.js)) but unused — a hook for a time axis.
- **Plot:** signal vs time; for k_obs-vs-[L] analysis, a secondary linear replot to extract
  k_on/k_off from the slope/intercept.
- **Current framework cost: medium.** New appmode, new `case` generating an exponential instead of
  a solved equilibrium, wire `axistype_time`. The fitter's grid search works on any `fun` returning
  `{d,m}`, so a time-course `fun` slots in, but the closures are hand-written.
- **Improved framework cost: low–medium.** `models/kinetics.js` with a `solve(p, t)` returning the
  exponential and `xAxis: "time"`. The fitter and renderer are model-agnostic. Main new piece is
  the linear-replot helper for k_obs analysis (a small model variant).

### 5.3 Inhibition: competitive / non-competitive / mixed
**Biophysics.** Michaelis–Menten with an inhibitor. The scaffolding is *already written but
stranded* at [calculations.js:821](../calculations.js): `alpha(inhib)`, `v0(S)`, `mm_rate_curve`,
`ic50_to_Ki`, `Ki_to_ic50`, and `enzyme.htm`. It needs finishing, not inventing.

- **Fix the model correctly:** the current `alpha()` only modifies `Km` (α) and treats
  non-/un-competitive incompletely. Proper forms:
  - competitive: v = Vmax·S / (α·Km + S), α = 1 + I/Ki
  - uncompetitive: v = Vmax·S / (Km + α′·S), α′ = 1 + I/Ki′
  - mixed: v = Vmax·S / (α·Km + α′·S)
  - non-competitive: mixed with Ki = Ki′
  (Vmax = kcat·E0.) The stub conflates these; the spec should implement all four with α and α′.
- **Parameters:** `kcat, Km, E0, I0, Ki, Ki′`, mechanism selector. Fit `kcat, Km, Ki` (and Ki′ for
  mixed). Cheng–Prusoff (`ic50_to_Ki`) already sketched for IC50↔Ki conversion (a README TODO).
- **Declare the globals:** `kcat, Km, E0, I0, Ki, inhib_mech` are currently undeclared — that is
  why the stub is inert.
- **Plot:** a family of v-vs-[S] curves at several [I]; optionally Lineweaver–Burk / Eadie–Hofstee
  linearisations as alternate x/y modes (mirrors the existing "alternative axis" pattern).
- **Current framework cost: medium.** Define `appmode_enzyme`, add sliders (the enzyme.htm stub
  lists them), declare globals, add an `update()` case drawing the curve family, add `fun_enzyme_*`
  to the fitter. The half-done files reduce the work.
- **Improved framework cost: low.** `models/inhibition.js` with `params` (incl. a discrete
  `mechanism` param) and `solve(p, S)` = the corrected `v0`. Curve *families* over [I] need the
  descriptor to express "one curve per value of parameter I" — a small generalisation of
  `buildCurves` worth designing once (also reused by §5.1's column families).

### 5.4 Two-state and three-state protein unfolding
**Biophysics.** Fraction folded vs denaturant. Two variants of the x-axis, both worth supporting:
chemical denaturation (linear extrapolation model, ΔG(D) = ΔG₀ − m·[D]) and thermal
(Gibbs–Helmholtz). Observed signal has **sloping native and unfolded baselines** — this is the
biophysically important detail the fit must include, not just the two-state Boltzmann.

- **Two-state:** K = exp(−ΔG/RT); f_U = K/(1+K). Observed = (y_N + s_N·[D]) + ((y_U + s_U·[D]) −
  (y_N + s_N·[D]))·f_U. Fit ΔG₀, m, and the four baseline params (y_N, s_N, y_U, s_U); Cm = ΔG₀/m.
- **Three-state (N ⇌ I ⇌ U):** two equilibria (ΔG₁, m₁, ΔG₂, m₂), populations f_N:f_I:f_U from
  the coupled equilibria, an intermediate signal level, and baselines. Populations here are exactly
  the "species populations" the project cares about — surface them in the data table/pie like the
  binding species.
- **Parameters:** ΔG(s), m-value(s), baselines. Fit all; report Cm/Tm.
- **Plot:** signal vs [denaturant] (or T); optionally the fraction-populated curves (f_N, f_I, f_U)
  as a second view — the population output the framework already knows how to render.
- **Current framework cost: medium.** New appmode(s); `update()` case computing populations and the
  baseline-modified signal; the y-axis is "signal" not "concentration", so the absolute/relative
  scaling logic needs a new mode. Two- and three-state can share a page with a selector.
- **Improved framework cost: low.** `models/unfolding.js` exporting two descriptors (or one with a
  `states: 2|3` param). `solve(p, D)` returns `{ d: [f_N, f_I, f_U], m: [1,1,1], signal }`; the
  populations flow into the existing data-table/pie path unchanged. Baselines are just more
  `params`.

### 5.5 Summary of add-cost

| model | current framework | improved framework | reuses existing |
|---|---|---|---|
| 5.1 competition matrix (plate) | **high** (no 2D concept) | medium | `calculate_ligands` cubic |
| 5.2 kinetics | medium | low–medium | `axistype_time` hook |
| 5.3 inhibition | medium | low | enzyme stub + `alpha/v0` |
| 5.4 unfolding 2-/3-state | medium | low | populations/pie path |

---

## 6. Order of work (recommended)

1. **Freeze current behaviour with golden tests** (§7) — capture curve arrays, populations, and fit
   results for all four models × all scale/x-axis/decadeshift combinations. Nothing else starts
   until this exists; it is the definition of "numerically identical".
2. **Extract pure cores without changing numerics:** move `expval`/solvers/math into modules,
   leave bodies byte-identical, keep the pages working. Re-run golden tests.
3. **Split `update()`** into `buildCurves` / `renderFigure` / `syncDataTable`; re-run golden tests.
4. **Introduce the model registry**; port the four models as descriptors wrapping the unchanged
   `calculate_*`. Re-run golden tests — this is the checkpoint that the refactor is behaviour-preserving.
5. **Centralise colours + add dark mode** (§3).
6. **Generalise the fitter** into `core/fit.js`; verify calcmode 1–4 outputs unchanged on the four
   models.
7. **Add new models** in the order 5.3 (cheapest, stub exists) → 5.4 → 5.2 → 5.1 (needs the plate
   renderer + global fit, the biggest new surface).

Each step is its own feature branch off `main` and is independently reversible.

---

## 7. Tests worth writing (noted, not yet written)

No test harness exists. Before writing tests, add a minimal runner (Vitest or plain `node --test`).
The suite falls into five groups.

**A. Golden / regression (the reproducibility guarantee — write these first).**
- For each of the four models, snapshot the full `curves[]` arrays for every combination of
  {absolute/relative/specificity y-scale} × {total/free/specificity x-mode} × a few `decadeshift`
  values, at fixed slider positions. These snapshots must not change across the entire refactor.
- Snapshot the data-table populations (conc, g/l, %) and pie values for fixed inputs.
- Snapshot `calculate_lsf` outputs (recovered slider positions / ext_values, and calcmode-4
  standard errors) for a fixed synthetic dataset per model.

**B. Solver / math unit tests.**
- Mass conservation: for every solver, Σ(species·multiplicity) equals the input totals
  (e.g. `calculate_ligand_total`: [PL]+[P] = E_0 and [PL]+[L] = S_0) across many random inputs.
- Limiting cases: K_D → 0 ⇒ fully bound; K_D → ∞ ⇒ unbound; half-saturation at [L]_free = K_D;
  homodimer high/low concentration limits.
- Numerical stability: the stable vs documented-unstable form of `calculate_ligand_total`
  ([calculations.js:49](../calculations.js)) agree in the well-conditioned regime and the stable
  one stays finite as [L] → 0.
- Cubic solvers: `solve_cubic_newton` and `solve_cubic_bisection` agree with each other and with a
  reference root for the competition polynomial; Newton's 20-cycle bail and bisection fallback are
  exercised; roots stay in [0, c_A].
- `calculate_ligands` competition: conservation of P, L, L′; recovers the single-ligand quadratic
  when the second ligand → 0.

**C. Slider / mapping (the numeric contract).**
- `expval` returns exact expected values at representative slider positions for `newdecadescale`
  true and false; monotonic in slider position; `decadeshift` shifts by exact powers of ten.
- Round-trip: value → nearest slider → value is stable; `valueinput_update` clamps out-of-range
  input as the current colour-coding implies.

**D. Fitter.**
- Recovers known parameters from noise-free synthetic data within tolerance for each model and each
  calcmode.
- calcmode-4 covariance: standard errors are positive, finite, and shrink as data density rises;
  Jacobian finite-difference matches an analytic derivative where one exists.
- Guard rails: "no data points", "no free parameters", "#points < #params", and the log-residual
  path all behave.
- `sum_squared_residuals` (linear and logarithmic) and `tinv_95` (table values n<20, approximation
  n≥20) against references.

**E. New-model tests (when built).**
- 5.1: per-well signal matches `calculate_ligands`; global fit recovers K_D, K_D′ from a synthetic
  plate; plate conservation per well.
- 5.2: exponential solve hits A_eq and t½ correctly; k_obs replot slope = k_on, intercept = k_off.
- 5.3: each mechanism reduces to plain MM at I=0; competitive raises apparent Km only; uncompetitive
  lowers both; Cheng–Prusoff IC50↔Ki round-trips.
- 5.4: two-state f_U = ½ at Cm; baseline-sloped signal reduces to plain Boltzmann when slopes = 0;
  three-state populations sum to 1 and reduce to two-state when the intermediate is destabilised.

Also worth a lightweight **DOM/render smoke test** (jsdom): `update()` runs without throwing for
each page and produces the expected number of `<polyline>`/legend nodes — cheap protection for the
render-layer split.

---

## 8. Branch & safety notes

- All work stays on feature branches off `main`; `main` is never the target. This planning doc is on
  `claude/binding-curve-refactor-5ad60d`.
- The refactor is explicitly behaviour-preserving for the four existing models; the golden suite
  (§7A) is the gate on every step.
- The biophysical read-outs the project values — per-species concentrations, g/l, %, pie
  populations, and the ΔG/K_A derived displays — are preserved and, for the new models (unfolding
  populations, plate signals), extended through the same populations/pie path rather than a parallel
  one.
