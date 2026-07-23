/*
	PROTEIN THERMODYNAMICS SIMULATIONS
	Copyright (C) 2021–2023 Johan Pääkkönen, Juha Rouvinen, University of Eastern Finland

	This program is free software; you can redistribute it and/or
	modify it under the terms of the GNU General Public License
	as published by the Free Software Foundation; either version 2
	of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

/*
	The GNU General Public License, version 2, is available at:
	(1) the file LICENSE in this folder
	(2) https://www.gnu.org/licenses/old-licenses/gpl-2.0.html
*/


// If you use this in your own research, please cite our paper: https://doi.org/10.1021/acsomega.2c00560

// binding-sim.js (refactor/06-embedding) — the <binding-sim> custom element.
//
// A drop-in embed: include this one classic <script> once, then place e.g.
//     <binding-sim model="ligand"></binding-sim>
// anywhere on a page. The element renders the corresponding simulation inside an
// isolated <iframe> pointing at the hosted static page. The iframe is a separate
// JS realm and CSS scope, so the sim's global state (S_0, curves, appmode, …) and
// the host page never collide — and multiple instances on one page are safe (each
// has its own globals). The sim pages ship as plain static files and are
// unchanged by this file.
//
// Attributes:
//   model   (required)  ligand | homodimer | ligands | receptors
//   ext                 boolean; loads the ?ext power-user features
//   theme               "dark" | "light" (best-effort, same-origin embeds only)
//   height              CSS length for the frame (default 900px; bare number → px)
//   width               CSS length (default 100%)
//   base                override the folder the sim pages are served from
//
// Isolation note: because each instance is an iframe, this is the safe embed even
// though the sims still use page globals internally. A future in-page Shadow-DOM
// component would require making that state per-instance first.

"use strict";

(function() {
	// Where the sim pages live: the folder containing this script (or data-base).
	var scriptEl = document.currentScript;
	var scriptBase = (scriptEl && scriptEl.dataset && scriptEl.dataset.base)
		? new URL(scriptEl.dataset.base, document.baseURI)
		: (scriptEl && scriptEl.src ? new URL(".", scriptEl.src) : new URL(".", document.baseURI));

	var MODELS = { ligand: 1, homodimer: 1, ligands: 1, receptors: 1, inhibition: 1 };

	function cssLen(v) {
		return /^[0-9.]+$/.test(v) ? v + "px" : v;
	}

	class BindingSim extends HTMLElement {
		static get observedAttributes() { return ["model", "ext", "theme", "height", "width", "base"]; }

		connectedCallback() {
			if(!this.shadowRoot) this.attachShadow({ mode: "open" });
			this.render();
		}

		attributeChangedCallback() {
			if(this.shadowRoot) this.render();
		}

		render() {
			var root = this.shadowRoot;
			var model = (this.getAttribute("model") || "").toLowerCase();
			var height = this.getAttribute("height") || "900px";
			var width = this.getAttribute("width") || "100%";
			var theme = this.getAttribute("theme");
			var base = this.getAttribute("base") ? new URL(this.getAttribute("base"), document.baseURI) : scriptBase;

			// Best-effort theme: the iframe (same origin) reads localStorage on load
			// via theme.js. Cross-origin embeds cannot be themed from the host — the
			// visitor can still use the toggle inside the embedded sim.
			if(theme === "dark" || theme === "light") {
				try { window.localStorage.setItem("bindingsims-theme", theme); } catch(e) {}
			}

			root.textContent = "";

			var style = document.createElement("style");
			style.textContent =
				":host { display: block; }" +
				".frame { width: " + cssLen(width) + "; height: " + cssLen(height) + "; }" +
				"iframe { width: 100%; height: 100%; border: 0; display: block; }" +
				".error { font: 14px/1.4 system-ui, sans-serif; color: #b00020; padding: 1em; }";
			root.appendChild(style);

			if(!MODELS[model]) {
				var err = document.createElement("div");
				err.className = "error";
				err.textContent = 'binding-sim: unknown model "' + model +
					'". Use one of: ligand, homodimer, ligands, receptors.';
				root.appendChild(err);
				return;
			}

			var url = new URL(model + ".htm", base);
			if(this.hasAttribute("ext")) url.search = "ext";

			var wrap = document.createElement("div");
			wrap.className = "frame";

			var iframe = document.createElement("iframe");
			iframe.src = url.href;
			iframe.setAttribute("title", model + " binding simulation");
			iframe.setAttribute("loading", "lazy");

			wrap.appendChild(iframe);
			root.appendChild(wrap);
		}
	}

	if(window.customElements && !window.customElements.get("binding-sim"))
		window.customElements.define("binding-sim", BindingSim);
})();
