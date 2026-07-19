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

// core/solvers.js — the binding math and cubic root solvers.
//
// Moved verbatim from calculations.js in refactor/01-core-math; the declarations
// below are byte-identical. `calculate_ligands_free` still reads `appmode` /
// `appmode_ligands` from the shared global scope (declared in simulation.js),
// exactly as it did when it lived in calculations.js — so behaviour is unchanged.

function calculate_ligand_free(E_0, S, K_D)
{
	var v = E_0 * S / (S + K_D);
	return { d: [v, E_0 - v], m: [1, 1] };
}

function calculate_ligand_total(E_0, S_0, K_D)
{
	/*
	// This the calculation described in the supporting information.
	// It is numerically unstable: when [L] -> 0 and [P] -> E_0,
	// precision is lost and [PL] becomes noisy.
	var b = (E_0 - S_0 + K_D) / 2;
	var L = -b + Math.sqrt(b * b + S_0 * K_D);
	var P = E_0 * K_D / (L + K_D);
	return { d: [E_0 - P, P, L], m: [1, 1, 0] };
	*/
	
	// This is the equivalent calculation but numerically stable.
	// When [PL] is calculated first, precision is never significantly lost
	// with possible parameter values in the simulation.
	var minus_b = (S_0 + E_0 + K_D) / 2;
	var v = minus_b - Math.sqrt(minus_b * minus_b - S_0 * E_0);
	return { d: [v, E_0 - v, S_0 - v], m: [1, 1, 0] };
}

function calculate_homodimer_free(P, K_D)
{
	return { d: [P * P / K_D, P], m: [2, 1] };
}

function calculate_homodimer(E_0, K_D)
{
	var P = (-K_D + Math.sqrt(K_D * K_D + 8 * K_D * E_0)) / 4;
	return { d: [P * P / K_D, P], m: [2, 1] };
}

function calculate_ligands_free(A, c_B, c_C, K_D, K_D2)
{
	var AB = c_B * A / (A + K_D);
	var AC = c_C * A / (A + K_D2);
	var B = c_B - AB;
	var C = c_C - AC;
	
	if(appmode === appmode_ligands)
		return { d: [AB, AC, A, B, C], m: [1, 1, 1, 0, 0] };
	else
		return { d: [AB, AC, B, C, A, AB / AC], m: [1, 1, 1, 1, 0, 0] };
}

function calculate_ligands(c_A, c_B, c_C, K_D, K_D2)
{
	var t1 = 1;
	var t2 = K_D2 + K_D + c_C + c_B - c_A
	var t3 = K_D2 * K_D + c_C * K_D + c_B * K_D2 - c_A * K_D2  - c_A * K_D;
	var t4 = -c_A * K_D2 * K_D;
	
	var A, B, C, AB, AC;
	
	A = solve_cubic_newton(t1, t2, t3, t4, c_A);
	
	if(!(A >= 0 && A < c_A))
		A = solve_cubic_bisection(t1, t2, t3, t4, 0, c_A);
	
	return calculate_ligands_free(A, c_B, c_C, K_D, K_D2);
}

function sum_squared_residuals(arr1, arr2, logarithmic)
{
	var i, r;
	var result = { residuals: [], sum: 0 };
	
	if(arr1.length !== arr2.length) return NaN;
	
	if(logarithmic)
	{
		for(i = 0; i < arr1.length; i++)
		{
			// r = Math.log(arr1[i].y) - Math.log(arr2[i].y);
			r = 1e-5 * Math.log(arr1[i].y / arr2[i].y);
			result.residuals.push(r);
			result.sum += r * r;
		}
	}
	else
	{
		for(i = 0; i < arr1.length; i++)
		{
			r = arr1[i].y - arr2[i].y;
			result.residuals.push(r);
			result.sum += r * r;
		}
	}
	
	return result;
}

function solve_cubic_newton(q1, q2, q3, q4, a0)
{
	// Cubic equation is solved using the Newton method.
	// This can fail to converge to the correct solution depending
	// on the choice of the initial guess.
	//
	// Equation: q1 * a^3 + q2 * a^2 + q3 * a + q4 = 0
	// initial guess: a0
	
	var a2, a1 = a0, cycles = 0;
	
	function val(x)
	{
		return (q1*x*x*x + q2*x*x + q3*x + q4);
		// return (q1*Math.pow(x,3) + q2*Math.pow(x,2) + q3*x + q4);
	}
	
	function der(x)
	{
		return (3*q1*x*x + 2*q2*x + q3);
		// return (3*q1*Math.pow(x,2) + 2*q2*x + q3);
	}
	
	// Abort if the derivative is very small because
	// the search would likely not converge
	
	if(val(a1) / der(a1) > 10) return NaN;
	
	while(1)
	{
		a2 = a1 - val(a1) / der(a1);
		
		if(Math.abs(a2 / a1 - 1) < 1e-3) return a2;
		a1 = a2;
		
		// Abort if too many cycles have been run
		// (the search will probably not converge)
		
		if(++cycles >= 20) return NaN; // TODO: try to find ways to reduce this
		
		// BTW: nice videos on the chaotic behaviour of the Newton method:
		// https://www.youtube.com/watch?v=-RdOwhmqP5s
		// https://www.youtube.com/watch?v=LqbZpur38nw
	}
}

function solve_cubic_bisection(q1, q2, q3, q4, amin, amax)
{
	// The root of a cubic equation is searched using the bisection method.
	//
	// Equation: q1 * a^3 + q2 * a^2 + q3 * a + q4 = 0
	// A root between amin and amax is searched.
	// Guaranteed to converge as long as a root exists.
	
	function val(x)
	{
		return (q1*x*x*x + q2*x*x + q3*x + q4);
		// return (q1*Math.pow(x,3) + q2*Math.pow(x,2) + q3*x + q4);
	}
	
	var xmin = amin;
	var xmax = amax;
	var ymin = val(xmin);
	var ymax = val(xmax);
	
	if(ymin === 0) return xmin;
	if(ymax === 0) return xmax;
	
	while(1)
	{
		var xmid = 0.5 * (xmin + xmax);
		var ymid = val(xmid);
		
		// This prevents an eternal loop in some edge cases
		if(xmid === xmin || xmid === xmax) return xmid;
		
		if(Math.abs(ymid) === 0)
		{
			return xmid;
		}
		else if(ymin * ymid > 0)
		{
			xmin = xmid;
			ymin = ymid;
		}
		else
		{
			xmax = xmid;
			ymax = ymid;
		}
		
		if(Math.abs(xmax / xmin - 1) < 1e-3) return xmid;
	}
}

function tinv_95(n)
{
	if(n < 0)
		return NaN;
	else if(n < 20)
		// first 19 accurate values from lookup table
		return ttable[n];
	else
		// my custom approximation
		// error zero at n=20, (n≈54.3) and n→∞
		//       positive when 20<n≤54
		//       negative when n≥55
		// maximum error about ±4.0×10⁻⁴ (±0.020%) at n=28 (+) and n=184 (-)
		return 1.95996 + 2.82680 * Math.pow(n, -1.03834);
}

var ttable = [
	Infinity, 12.70620, 4.30265, 3.18245, 2.77645,
	2.57058,   2.44691, 2.36462, 2.30600, 2.26216,
	2.22814,   2.20099, 2.17881, 2.16037, 2.14479,
	2.13145,   2.11991, 2.10982, 2.10092, 2.09302
];

export {
	calculate_ligand_free, calculate_ligand_total,
	calculate_homodimer_free, calculate_homodimer,
	calculate_ligands_free, calculate_ligands,
	sum_squared_residuals, solve_cubic_newton, solve_cubic_bisection, tinv_95,
};
