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

// theme.js (refactor/05) — light/dark toggle for the page chrome.
//
// Loaded as a CLASSIC script in <head> so the saved theme is applied before the
// first paint (no flash of the wrong theme — deferred modules would flash). Only
// the surrounding page is themed, via the CSS custom properties in style.css;
// the plot figure stays a light "card" in both themes, so the exported SVG and
// the golden tests are unaffected. Default is light: the white look is preserved
// unless the visitor opts into dark.

"use strict";

(function() {
	var KEY = "bindingsims-theme";

	function apply(theme) {
		if(theme === "dark")
			document.documentElement.setAttribute("data-theme", "dark");
		else
			document.documentElement.removeAttribute("data-theme");
	}

	function current() {
		return document.documentElement.getAttribute("data-theme") === "dark" ? "dark" : "light";
	}

	function update_button() {
		var btn = document.getElementById("theme_toggle");
		if(!btn) return;
		var dark = current() === "dark";
		btn.textContent = dark ? "☀ Light mode" : "☾ Dark mode";
		btn.setAttribute("aria-pressed", dark ? "true" : "false");
	}

	// Apply the saved theme immediately (before paint).
	var saved = null;
	try { saved = localStorage.getItem(KEY); } catch(e) {}
	apply(saved);

	// Exposed for the inline onclick on the toggle button.
	window.toggle_theme = function() {
		var next = current() === "dark" ? "light" : "dark";
		apply(next);
		try { localStorage.setItem(KEY, next); } catch(e) {}
		update_button();
	};

	if(document.readyState === "loading")
		document.addEventListener("DOMContentLoaded", update_button);
	else
		update_button();
})();
