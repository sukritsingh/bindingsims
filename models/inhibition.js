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

// models/inhibition.js (feature/07-inhibition) — reversible enzyme inhibition.
//
// Steady-state Michaelis–Menten with an inhibitor. The general mixed form is
//
//     v = Vmax·[S] / (α·Km + α′·[S]),   Vmax = kcat·[E]0
//     α  = 1 + [I]/Ki      (inhibitor binding to free enzyme  E)
//     α′ = 1 + [I]/Ki′     (inhibitor binding to the ES complex)
//
// and the four textbook mechanisms are special cases:
//   competitive     α′ = 1           (Ki′ → ∞; raises apparent Km, Vmax unchanged)
//   uncompetitive   α  = 1           (Ki  → ∞; lowers apparent Km and Vmax together)
//   non-competitive Ki′ = Ki → α = α′ (lowers Vmax, apparent Km unchanged)
//   mixed           independent Ki, Ki′
//
// This replaces the stranded stub in calculations.js, whose alpha() only ever
// modified Km and so conflated the uncompetitive / non-competitive / mixed cases.

/**
 * @typedef {Object} InhibParams
 * @property {number} kcat        turnover number (s⁻¹)
 * @property {number} E0          total enzyme concentration (mol l⁻¹)
 * @property {number} Km          Michaelis constant (mol l⁻¹)
 * @property {number} Ki          inhibition constant for free enzyme (mol l⁻¹)
 * @property {number} KiP         inhibition constant for the ES complex (mol l⁻¹)
 * @property {string} mechanism   "none" | "competitive" | "uncompetitive" | "noncompetitive" | "mixed"
 */

/** α, α′ multipliers for a given inhibitor concentration and mechanism. */
function factors(p, I) {
	switch(p.mechanism) {
		case "competitive":    return { a: 1 + I / p.Ki,  aP: 1 };
		case "uncompetitive":  return { a: 1,             aP: 1 + I / p.KiP };
		case "noncompetitive": return { a: 1 + I / p.Ki,  aP: 1 + I / p.Ki };
		case "mixed":          return { a: 1 + I / p.Ki,  aP: 1 + I / p.KiP };
		default:               return { a: 1,             aP: 1 };   // "none": no inhibitor
	}
}

/**
 * Initial velocity v at substrate [S] and inhibitor [I].
 * @param {InhibParams} p
 * @param {number} S   substrate concentration (mol l⁻¹)
 * @param {number} I   inhibitor concentration (mol l⁻¹)
 * @returns {number}   velocity (mol l⁻¹ s⁻¹)
 */
export function enzyme_rate(p, S, I) {
	var f = factors(p, I);
	var Vmax = p.kcat * p.E0;
	return Vmax * S / (f.a * p.Km + f.aP * S);
}

/**
 * IC50 — the inhibitor concentration giving half the uninhibited velocity at [S]
 * (Cheng–Prusoff). General form: IC50 = ([S] + Km) / (Km/Ki + [S]/Ki′), which
 * reduces to the familiar per-mechanism expressions.
 * @param {InhibParams} p
 * @param {number} S
 * @returns {number}
 */
export function ic50(p, S) {
	switch(p.mechanism) {
		case "competitive":    return p.Ki * (1 + S / p.Km);
		case "uncompetitive":  return p.KiP * (1 + p.Km / S);
		case "noncompetitive": return p.Ki;
		case "mixed":          return (S + p.Km) / (p.Km / p.Ki + S / p.KiP);
		default:               return Infinity;   // "none": never inhibited
	}
}

/** Ki from a measured IC50 for a competitive inhibitor (classic Cheng–Prusoff). */
export function ic50_to_Ki(ic50val, S, Km) {
	return ic50val / (1 + S / Km);
}

/** IC50 expected for a competitive inhibitor of a given Ki. */
export function Ki_to_ic50(Ki, S, Km) {
	return Ki * (1 + S / Km);
}

// Slider index → parameter, shared by the page wiring and the fitter dispatch:
//   2 Km, 4 kcat, 5 E0, 3 S (substrate, swept), 11 I0, 12 Ki, 13 Ki′.
// The mechanism is the ambient `inhib_mech` global (like appmode for the cubic
// solver): read here so the fitter and the figure agree on the active form.
function paramsFromM(m) {
	return {
		kcat: m[4], E0: m[5], Km: m[2], Ki: m[12], KiP: m[13],
		mechanism: (typeof inhib_mech !== "undefined") ? inhib_mech : "none",
	};
}

export const inhibition = {
	id: "inhibition",
	legacyAppmode: 4,          // appmode_enzyme (see simulation.js)
	datalabels: [],            // a rate model: no species populations / pie / table
	axisLabels(a) {
		var { xscale_alternative, scale_absolute } = a;
		return {
			x: "Substrate concentration (mol l\uEEE1\u22121\uEEE0)",
			y: scale_absolute
				? "Reaction rate (mol l\uEEE1\u22121\uEEE0 s\uEEE1\u22121\uEEE0)"
				: "Relative reaction rate (v / V\uEEE2max\uEEE0)",
		};
	},
	// Fit v-vs-[S] data at the current inhibitor concentration [I]0 (m[11]).
	fitSolve(m, x, xscale_alternative) {
		var p = paramsFromM(m);
		return { d: [enzyme_rate(p, x, m[11])], m: [1] };
	},
};
