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

// core/fit.js — the DOM-free numerical core of the least-squares fitter.
//
// Extracted from calculate_lsf() (refactor/04). calculate_lsf stays in
// calculations.js as the thin DOM adapter: it reads the sliders / checkboxes /
// radios into a plain `cfg` object, calls these functions, and writes the
// results back to the DOM. The arithmetic (grid search + calcmode-4 Jacobian /
// Gauss-Jordan covariance) is unchanged and is model-agnostic: the model enters
// only through `cfg.predict(m, x) -> { d, m }` (the descriptor's fitSolve bound
// to the current x-axis mode). Verified by the fit goldens (calcmodes 1–4).

import { expval } from './expval.js';
import { sum_squared_residuals } from './solvers.js';

/**
 * @typedef {Object} FitConfig
 * @property {{x:number,y:number}[]} data           the fitted data points
 * @property {(m:number[], x:number)=>{d:number[],m:number[]}} predict  model solve
 * @property {number[]} m                            full parameter vector (by slider index)
 * @property {number[]} sliderIndices                slider index of each free parameter
 * @property {{min:number,max:number,value:number,ext_value:number}[]} bounds  per free param
 * @property {number} calcmode                       1..4
 * @property {number} co                             observed-species column
 * @property {number} scale_absolute                 0 relative, 1 absolute, 2 specificity
 * @property {number} appmode
 * @property {number} appmode_receptors
 * @property {boolean} logresid
 * @property {Array} expparams
 */

/**
 * Coarse-to-fine grid search over the free parameters.
 * @param {FitConfig} cfg
 * @returns {{kp:number[], resid_min:number[], srsum_min:number}}
 */
export function gridSearch(cfg) {
	var data = cfg.data;
	var predict = cfg.predict;
	var m = cfg.m;
	var slider_indices = cfg.sliderIndices;
	var bounds = cfg.bounds;
	var calcmode = cfg.calcmode;
	var co = cfg.co;
	var scale_absolute = cfg.scale_absolute;
	var appmode = cfg.appmode;
	var appmode_receptors = cfg.appmode_receptors;
	var logresid = cfg.logresid;
	var expparams = cfg.expparams;

	var n_free = bounds.length;
	var kmin = [], kmax = [], kstep = [], kp = [], k = [];
	var theor = [];
	var srsum_min, resid_min, fitdata;
	var i, j;

	for(j = 0; j < (calcmode === 1 ? 2 : 1); j++)
	{
		if(j === 0)
		{
			for(i = 0; i < n_free; i++)
			{
				switch(calcmode)
				{
					case 1:
						kmin[i] = bounds[i].min;
						kmax[i] = bounds[i].max;
						kstep[i] = 10;
						break;
					case 2:
						kmin[i] = bounds[i].min;
						kmax[i] = bounds[i].max;
						kstep[i] = 1;
						break;
					case 3:
						kmin[i] = Math.max(bounds[i].min, bounds[i].value - 5);
						kmax[i] = Math.min(bounds[i].max, bounds[i].value + 5);
						kstep[i] = 1;
						break;

					// extmode only
					case 4:
						kstep[i] = bounds[i].ext_value * 0.0001;
						kmin[i] = bounds[i].ext_value - kstep[i] * 10;
						kmax[i] = bounds[i].ext_value + kstep[i] * 10.1;
						break;
				}
			}
		}
		else
		{
			// calcmode is always 1 in this case

			for(i = 0; i < n_free; i++)
			{
				kmin[i] = Math.max(bounds[i].min, kp[i] - 20);
				kmax[i] = Math.min(bounds[i].max, kp[i] + 20);
				kstep[i] = 1;
			}
		}

		for(i = 0; i < n_free; i++)
		{
			k[i] = kmin[i];
		}

		var loop = true;

		srsum_min = Number.MAX_VALUE;
		resid_min = undefined;

		while(loop)
		{
			for(i = 0; i < n_free; i++)
			{
				if(calcmode === 4)
				{
					m[slider_indices[i]] = k[i];
				}
				else
				{
					m[slider_indices[i]] = expval(k[i],
						expparams[slider_indices[i]][0],
						expparams[slider_indices[i]][1],
						expparams[slider_indices[i]][2]);
				}
			}

			for(i = 0; i < data.length; i++)
			{
				var cd = predict(m, data[i].x);
				var divisor;

				if(scale_absolute || appmode === appmode_receptors)
					divisor = 1;
				else
					divisor = cd.d.reduce(function(t,v,ii){ return t + v * cd.m[ii]; }, 0);

				theor[i] = {
					x: data[i].x,
					y: cd.d[co] / divisor * (scale_absolute || appmode === appmode_receptors ? 1 : cd.m[co])
				};
			}

			fitdata = sum_squared_residuals(data, theor, logresid);

			if(fitdata.sum < srsum_min)
			{
				srsum_min = fitdata.sum;
				resid_min = fitdata.residuals;
				kp = k.slice();
			}

			for(i = 0; i < n_free; i++)
			{
				if((k[i] += kstep[i]) <= kmax[i])
					break;
				else if(i < n_free - 1)
					k[i] = kmin[i];
				else
					loop = false;
			}
		}
	}

	// kstep is returned because the calcmode-4 convergence test in the DOM adapter
	// compares |kp - ext_value| against kstep * 0.1.
	return { kp: kp, resid_min: resid_min, srsum_min: srsum_min, kstep: kstep };
}

/**
 * Calcmode-4 variance–covariance matrix: finite-difference Jacobian, then
 * (J'J)^-1 * mse via Gauss-Jordan elimination.
 * @param {FitConfig & {srsum_min:number, resid_min:number[]}} cfg
 * @param {number[]} kp   best-fit parameter values (calcmode-4 units)
 * @returns {number[]|undefined}  row-major J_cols×J_cols covariance, or undefined
 */
export function covariance(cfg, kp) {
	var data = cfg.data;
	var predict = cfg.predict;
	var m = cfg.m;
	var slider_indices = cfg.sliderIndices;
	var bounds = cfg.bounds;
	var co = cfg.co;
	var scale_absolute = cfg.scale_absolute;
	var appmode = cfg.appmode;
	var appmode_receptors = cfg.appmode_receptors;
	var logresid = cfg.logresid;
	var srsum_min = cfg.srsum_min;
	var resid_min = cfg.resid_min;

	var n_free = bounds.length;
	var covb;
	var i, j, n;

	if(resid_min.length - n_free > 0)
	{
		// calculate mean square error (mse)

		let mse = srsum_min / (resid_min.length - n_free);

		// calculate jacobian matrix J

		let J_rows = resid_min.length;
		let J_cols = n_free;
		let J = new Array(J_rows * J_cols);

		for(i = 0; i < n_free; i++)
		{
			m[slider_indices[i]] = kp[i];
		}

		for(i = 0; i < data.length; i++)
		{
			for(j = 0; j < n_free; j++)
			{
				let delta = bounds[j].ext_value * 0.0001;
				let temp = m[slider_indices[j]];

				m[slider_indices[j]] = temp - delta;
				var cd1 = predict(m, data[i].x);

				m[slider_indices[j]] = temp + delta;
				var cd2 = predict(m, data[i].x);

				m[slider_indices[j]] = temp;

				var divisor1, divisor2;

				if(scale_absolute || appmode === appmode_receptors)
				{
					divisor1 = 1;
					divisor2 = 1;
				}
				else
				{
					divisor1 = cd1.d.reduce(function(t,v,ii){ return t + v * cd1.m[ii]; }, 0);
					divisor2 = cd2.d.reduce(function(t,v,ii){ return t + v * cd2.m[ii]; }, 0);
				}

				let y1 = cd1.d[co] / divisor1 * (scale_absolute || appmode === appmode_receptors ? 1 : cd1.m[co]);
				let y2 = cd2.d[co] / divisor2 * (scale_absolute || appmode === appmode_receptors ? 1 : cd2.m[co]);

				if(logresid)
					J[i * J_cols + j] = 1e-5 * Math.log(y2 / y1) / (2 * delta);
				else
					J[i * J_cols + j] = (y2 - y1) / (2 * delta);
			}
		}

		// calculate J'*J, make an augmented matrix for Gauss-Jordan

		let JJ = new Array(J_cols * J_cols * 2);

		for(j = 0; j < J_cols; j++)
		for(i = 0; i < J_cols; i++)
		{
			let t = 0;

			for(n = 0; n < J_rows; n++)
			{
				t += J[n * J_cols + i] * J[n * J_cols + j];
			}

			JJ[j * 2 * J_cols + i] = t;
			JJ[j * 2 * J_cols + i + J_cols] = (i == j) ? 1 : 0;
		}

		// invert J'*J using Gauss-Jordan elimination
		// (this is surprisingly simple, though it could be written to be clearer)
		// (this may be numerically unstable, but I have not yet run into major issues)

		function subtract_row(r1, r2, f)
		{
			let i;
			for(i = 0; i < 2 * J_cols; i++)
			{
				JJ[r1 * 2 * J_cols + i] -= f * JJ[r2 * 2 * J_cols + i];
			}
		}

		for(i = 0; i < J_cols; i++)
		{
			for(j = i; j < J_cols; j++)
			{
				subtract_row(j, i, (JJ[j * 2 * J_cols + i] - (i == j ? 1 : 0)) / JJ[i * 2 * J_cols + i]);
			}
		}

		for(i = J_cols - 1; i >= 0; i--)
		{
			for(j = i - 1; j >= 0; j--)
			{
				subtract_row(j, i, JJ[j * 2 * J_cols + i] / JJ[i * 2 * J_cols + i]);
			}
		}

		// calculate the variance-covariance matrix

		covb = new Array(J_cols * J_cols);

		for(j = 0; j < J_cols; j++)
		for(i = 0; i < J_cols; i++)
		{
			covb[j * J_cols + i] = JJ[j * 2 * J_cols + i + J_cols] * mse;
		}
	}

	return covb;
}
