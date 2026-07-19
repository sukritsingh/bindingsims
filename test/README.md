# Tests — the golden gate

This suite exists to guarantee the refactor stays **numerically identical**. It drives the
**unmodified** simulation code inside jsdom (a headless browser), captures its outputs, and asserts
they never change. See [../docs/BLUEPRINTS.md](../docs/BLUEPRINTS.md) §5–§6 for the plan.

## Run

```bash
npm install      # dev-only: jsdom + vitest. Nothing here ships to users.
npm test         # run the suite once
npm run test:watch
```

## What is checked

| File | Locks |
|---|---|
| `golden.test.js` → `golden/math.golden.json` | every solver's arithmetic + the `expval` slider mapping, over fixed input grids |
| `golden.test.js` → `golden/curves.golden.json` | the full figure data (curves, pie, axes, data table) for all 4 models × their scenarios |
| `golden.test.js` → `golden/fit.golden.json` | `calculate_lsf` recovered slider positions from noise-free synthetic data |
| `invariants.test.js` | analytic truths independent of the snapshots: mass conservation, half-saturation, cubic-solver agreement, K_D recovery |

## How it works

- `helpers/loadSim.js` — loads a page into jsdom, injects `simulation.js`/`calculations.js`/`update.js`
  as inline scripts (so their globals populate `window`, exactly as in a browser), and runs
  `preinit()`/`init()`. **This is the only file that changes when the refactor moves code into
  modules** — the golden JSON stays frozen.
- `helpers/scenarios.js` — scenario definitions + `snapshot()` + the deterministic serializers.
- `helpers/mathCases.js` — input grids for the solvers and the fitter cases.
- `helpers/compare.js` — tolerant structural comparison. Numbers are matched within a
  relative tolerance (`1e-3`) rather than bit-exactly. Two things force this across
  platforms (macOS dev vs Linux CI): `Math.pow`/`log10` differ by ~1 ULP, and a few
  curve points in the free-protein tail (`[P] = E_0 - v`) are catastrophic-cancellation
  garbage that isn't bit-reproducible (observed delta ~1.7e-4 on ~11 isolated points).
  A real behaviour change shifts values by whole percent — far above the tolerance — so
  it's still caught; the analytic invariants stay tight (`1e-12`) as a second check.

## Golden discipline

The golden JSON is the frozen reference. **Never edit it by hand, and do not regenerate it during a
refactor step** — a diff there means behaviour changed, which is exactly what the gate should catch.
Only when you *intend* to change behaviour (e.g. adding a new model in a feature branch) do you run:

```bash
npm run golden   # re-captures all three golden files from the current code
```

and review the diff before committing.
