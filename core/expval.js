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

// core/expval.js — the slider<->value contract (the frozen numerics).
//
// Moved verbatim from simulation.js in refactor/01-core-math; the declarations
// below are byte-identical. `newdecadescale` and `decadeshift` are still read
// from the shared global scope (declared in simulation.js), exactly as they
// resolved when these functions lived there — so behaviour is unchanged.

var decadetable = [
	10, 11, 12, 13, 14, 15, 16, 17, 18, 20,
	22, 24, 26, 28, 30, 32, 34, 36, 38, 40,
	45, 50, 55, 60, 65, 70, 75, 80, 85, 90
];

function expval(v, e, minexp, maxexp)
{
	if(newdecadescale)
		return decadetable[v % 30] * Math.pow(10, Math.floor(v / 30) - 10 + decadeshift);
	else
		return Math.pow(e, minexp + v * (maxexp - minexp) / 240 + decadeshift);
}

var expparams = [
	/*  0 */ [0, 0, 0],
	/*  1 */ [10, -1, 5],
	/*  2 */ [10, -7, -1],
	/*  3 */ [10, -9, -1],
	/*  4 */ [10, 0, 6],
	/*  5 */ [10, -9, -1],
	/*  6 */ [10, -5, -1],
	/*  7 */ [10, -9, -1],
	/*  8 */ [0, 0, 0],
	/*  9 */ [10, -9, -1],
	/* 10 */ [10, -9, -1],
	/* 11 */ [10, -9, 1],
	/* 12 */ [10, -9, -1],
	/* 13 */ [10, -9, -1],
	/* 14 */ [10, -12, -4],
	/* 15 */ [10, -12, -4]
];

function expval_wrap(v, index)
{
	return expval(v, expparams[index][0], expparams[index][1], expparams[index][2]);
}

export { decadetable, expval, expparams, expval_wrap };
