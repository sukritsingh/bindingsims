# Enzyme inhibition model

The enzyme-inhibition simulation ([`inhibition.htm`](../inhibition.htm)) is the
first feature model (branch `feature/07-inhibition`). It finishes the enzyme stub
that was stranded in `calculations.js`, whose `alpha()` only ever modified `Km`
and so conflated the mechanisms.

## Biophysics

Steady-state Michaelis–Menten with a reversible inhibitor. The general **mixed**
form is

```
v = Vmax·[S] / (α·Km + α′·[S]),   Vmax = kcat·[E]0
α  = 1 + [I]/Ki      (inhibitor binding free enzyme  E)
α′ = 1 + [I]/Ki′     (inhibitor binding the ES complex)
```

The four textbook mechanisms are special cases:

| Mechanism | α | α′ | Effect |
|---|---|---|---|
| competitive | 1 + [I]/Ki | 1 | raises apparent Km, Vmax unchanged |
| uncompetitive | 1 | 1 + [I]/Ki′ | lowers apparent Km and Vmax together |
| non-competitive | 1 + [I]/Ki | = α (Ki′ = Ki) | lowers Vmax, Km unchanged |
| mixed | 1 + [I]/Ki | 1 + [I]/Ki′ | independent Ki, Ki′ |

**Cheng–Prusoff.** The IC50 (the [I] giving half the uninhibited velocity at the
current [S]) is `IC50 = ([S] + Km) / (Km/Ki + [S]/Ki′)`, which reduces to the
familiar per-mechanism forms (e.g. competitive `Ki·(1 + [S]/Km)`). It is shown
live on the page and is exposed as `ic50()` / `ic50_to_Ki()` / `Ki_to_ic50()`.

The math and these limits are pinned in [`test/inhibition.test.js`](../test/inhibition.test.js).

## Parameters and controls

Sliders: `[E]0` (enzyme), `[S]` (substrate reference/marker), `kcat`, `Km`,
`[I]0` (inhibitor), `Ki`, `Ki′`, plus a mechanism selector. The plot is a family
of v-vs-[S] curves at several `[I]` (0, ½, 1, 2, 4 × `[I]0`); the vertical scale
toggles absolute rate vs relative (v / Vmax) and the horizontal axis toggles
log/linear [S]. The fitter recovers `kcat`, `Km`, `Ki` (and `Ki′`) from a
v-vs-[S] dataset via the model registry (`models/inhibition.js` `fitSolve`).

## How it is wired (current, imperative framework)

This model was added the same way as the four binding models rather than as a
pure descriptor, because the generic figure builder / `family` layout were
deferred during the refactor:

- `models/inhibition.js` — the corrected rate math, IC50, and the registry
  descriptor (`legacyAppmode: 4`, empty `datalabels`, `axisLabels`, `fitSolve`).
- `simulation.js` — `appmode_enzyme`, the `kcat/Km/I0/Ki/Ki2/inhib_mech` globals,
  `slider_input` cases (2 Km, 4 kcat, 11 I0, 12 Ki, 13 Ki′; [S]/[E]0 reuse 3/5),
  and `mechanism_input()`.
- `update.js` — a `buildFigure` case that draws the v-vs-[S] family (no pie / no
  species table); `renderFigure` is reused unchanged.
- `core/bootstrap.js` / `test/helpers/loadSim.js` expose `enzyme_rate` / `ic50`
  and (for tests) bridge the `inhib_mech` global.

It is embeddable like the others: `<binding-sim model="inhibition"></binding-sim>`.

## Known limitation

The concentration sliders use the app's decade scale (~10⁻¹⁰ – 10⁰), so **kcat is
capped near ~1 s⁻¹** rather than the realistic 1–10⁶ s⁻¹. Because `Vmax = kcat·[E]0`
only scales the vertical axis, the inhibition *shapes and patterns* (the teaching
content) and the *relative-rate* view are unaffected; only the absolute rate
magnitude is small. A dedicated rate-scale slider is a possible follow-up.
