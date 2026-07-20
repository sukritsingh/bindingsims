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

// core model descriptor: receptors

import { calculate_ligands, calculate_ligands_free } from '../core/solvers.js';

export const receptors = {
	id: "receptors",
	legacyAppmode: 1,
	datalabels: ["[PL]", "[P\u2032L]", "[P]", "[P\u2032]"],
	axisLabels(a) {
		var { xscale_alternative, scale_absolute, xmagnitude, ymagnitude, magnitude_string } = a;
		var xlabelstr = "", ylabelstr = "";
			xlabelstr = (xscale_alternative ? "Free" : "Total") + " ligand L concentration (" + magnitude_string(xmagnitude) + "mol l\uEEE1\u22121\uEEE0)";
			// ylabelstr = scale_absolute ? "Protein species concentration (mol l\uEEE1\u22121\uEEE0)" : "Protein species concentration (" + magnitude_string(ymagnitude) + "mol l\uEEE1\u22121\uEEE0)";
			switch(scale_absolute)
			{
				case 2: ylabelstr = "Receptor P specificity" + (ymagnitude ? " / " + magnitude_string(ymagnitude, true) : ""); break;
				case 1: ylabelstr = "Protein species concentration (mol l\uEEE1\u22121\uEEE0)"; break;
				case 0: ylabelstr = "Protein species concentration (" + magnitude_string(ymagnitude) + "mol l\uEEE1\u22121\uEEE0)"; break;
			}
		return { x: xlabelstr, y: ylabelstr };
	},
	fitSolve(m, x, xscale_alternative) {
		return xscale_alternative
			? calculate_ligands_free(x, m[10], m[3], m[7], m[9])
			: calculate_ligands(x, m[10], m[3], m[7], m[9]);
	},
};
