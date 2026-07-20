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

// core/models registry — the single lookup for every model descriptor.
// New models register here; the model-contract test (test/models.test.js) runs
// generic invariants over everything exported from this file.

import { ligand } from './ligand.js';
import { homodimer } from './homodimer.js';
import { ligands } from './competing-ligands.js';
import { receptors } from './competing-receptors.js';

export const models = { ligand, homodimer, ligands, receptors };

const byAppmode = {};
for (const m of Object.values(models)) byAppmode[m.legacyAppmode] = m;

/** Look up a descriptor by the legacy integer appmode id (see simulation.js). */
export function modelByAppmode(appmode) { return byAppmode[appmode]; }

// datalabels indexed by legacy appmode id, mirroring the old global that used to
// live in calculations.js — now derived from the registry (single source).
export const datalabels = [];
for (const m of Object.values(models)) datalabels[m.legacyAppmode] = m.datalabels;
