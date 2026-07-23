# Embedding bindingsims

The simulations are plain static files, so embedding them in another site (a lab
page, a course, a blog post) is easy. There are three ways to do it, cheapest to
cleanest. All of them require the files to be **served over HTTP(S)** — the sims
use native ES modules, which browsers refuse to load from a `file://` path. Any
real host (GitHub Pages, a lab server, an intranet) already satisfies this; for
local testing run `python3 -m http.server` in the repo and use the
`http://localhost:…` URLs.

## 1. `<iframe>` — zero code

Point an iframe at a hosted simulation page:

```html
<iframe src="https://your-host/bindingsims/ligand.htm"
        width="100%" height="900" title="Ligand binding simulation"></iframe>
```

Fully isolated (its own JS and CSS), works everywhere, no scripts to include.
Downside: you size it by hand and it reads as a bolted-on box. Pick the page for
the model you want — `ligand.htm`, `homodimer.htm`, `ligands.htm`,
`receptors.htm` — and append `?ext` for the power-user features. See
[`examples/iframe.html`](../examples/iframe.html).

## 2. `<binding-sim>` web component — one script, then a tag

Include the component once, then drop the tag anywhere:

```html
<script src="https://your-host/bindingsims/binding-sim.js"></script>

<binding-sim model="ligand"></binding-sim>
<binding-sim model="receptors" ext height="1000"></binding-sim>
```

The element renders the sim inside an isolated iframe, so **multiple instances on
one page are safe** and nothing collides with the host page. It resolves the sim
pages relative to `binding-sim.js`'s own URL, so host that script alongside the
`.htm` files (or set `data-base` on the script / `base` on the element to point
elsewhere).

| Attribute | Values | Notes |
|---|---|---|
| `model` | `ligand` \| `homodimer` \| `ligands` \| `receptors` | required |
| `ext` | boolean (presence) | loads the `?ext` power-user features |
| `theme` | `dark` \| `light` | best-effort, **same-origin embeds only** (sets shared `localStorage`); otherwise use the toggle inside the sim |
| `height` | CSS length | frame height, default `900px`; a bare number means pixels |
| `width` | CSS length | default `100%` |
| `base` | URL | override where the sim pages are served from |

See [`examples/web-component.html`](../examples/web-component.html).

### Why an iframe and not in-page Shadow DOM?

The tag is deliberately backed by an iframe. The sims still keep their state in
page globals (`S_0`, `curves`, `appmode`, …) and query the DOM by fixed ids — a
pragmatic choice from the refactor. An iframe gives each instance its own JS
realm, which is what makes multiple safe instances possible today. A future
in-page Shadow-DOM version would first require making that state per-instance and
rescoping every DOM lookup.

## 3. Folder drop-in

Because every path in the sims is relative, you can copy a sim's files
(`*.htm`, `*.js`, `core/`, `models/`, `style.css`, `theme.js`) into a subfolder
of another **static** site and link to the page directly. Serve over HTTP, as
above. This is the lowest-ceremony way to self-host without an iframe, and is
safe now that the shared code lives in ES modules under `core/` and `models/`.

## Themes

Each sim has its own light/dark toggle (top-right) and remembers the choice in
`localStorage` per origin. Same-origin embeds inherit that choice automatically;
the `theme` attribute above is a convenience for setting it from the host.
