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

// core/bootstrap.js — the compatibility shim (refactor/01-core-math, decision:
// "window shim"). The page still loads simulation.js / calculations.js /
// update.js as classic scripts, so their functions and state stay on `window`
// exactly as before. This module imports the core ES modules and re-exposes
// their functions on `window` under the same names, so the classic code can go
// on calling `expval(...)`, `calculate_ligands(...)`, `render_text_svg(...)`
// etc. by bare name. The moved functions in turn read the shared state
// (`newdecadescale`, `decadeshift`, `svg_xmlns`, `delete_all_children`,
// `appmode`) off `window`, precisely as they did when they lived in those
// files — so nothing observable changes. A later step retires these globals as
// the page itself becomes modular.
//
// Module scripts are deferred, so this runs after the classic scripts have
// defined their globals and before the page's `load` event fires `init()`.

import * as expvalMod from './expval.js';
import * as solversMod from './solvers.js';
import * as formatMod from './format.js';

Object.assign(window, expvalMod, solversMod, formatMod);
