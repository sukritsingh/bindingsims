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

"use strict";

// `datalabels` and the per-model solver dispatch now come from the model
// registry (models/index.js), exposed on window by core/bootstrap.js. See the
// use of modelByAppmode() in calculate_lsf below.

function calculate_lsf()
{
	if(datapoints.length < 1)
	{
		var label = document.getElementById("calcstatus");
		if(label)
		{
			label.innerHTML = "There are no data points.";
			label.style.color = "red";
		}
		
		return;
	}
	
	var k1, k2, k3, k4;
	var m1, m2, m3;
	var k1p, k2p, k3p, k4p;
	var m1p, m2p, m3p;
	
	var k1min, k1max, k1step;
	var k2min, k2max, k2step;
	var k3min, k3max, k3step;
	var k4min, k4max, k4step;
	
	var srsum_min;
	var resid_min;
	var i, j, n;
	var co;
	
	var calcmode;
	
	var logresid = extmode && document.getElementById("ext_calcmode_log").checked;
	
	for(i = 1; i <= 4; i++)
	{
		var ele = document.getElementById("calcmode" + i);
		
		if(ele && ele.checked)
		{
			calcmode = i;
			break;
		}
	}
	
	var slider_indices = [];
	var sliders = [];
	var kp;
	var m = [];

	// Solver dispatch comes from the model registry (models/*.js). predict(m, x)
	// evaluates the model's fitSolve at a data point; the numerical fit itself
	// (grid search + covariance) lives in core/fit.js and never touches the DOM.
	var model = modelByAppmode(appmode);
	var predict = function(mm, x) { return model.fitSolve(mm, x, xscale_alternative); };
	
	var fixvallist = document.getElementsByClassName("fixval");
	
	for(i = 0; i < fixvallist.length; i++)
	{
		var index = Number(fixvallist[i].id.substring(6));
		var s = document.getElementById("slider" + index);
		
		m[index] = expval(Number(s.value), expparams[index][0], expparams[index][1], expparams[index][2]);
		
		if(!(fixvallist[i].disabled || !fixvallist[i].checked))
		{
			slider_indices.push(index);
			sliders.push(s);
			
			// only in extmode
			if(calcmode === 4 && !(s.ext_value))
			{
				s.ext_value = m[index];
			}
		}
	}
	
	if(sliders.length < 1)
	{
		var label = document.getElementById("calcstatus");
		if(label)
		{
			label.innerHTML = "There are no free parameters to optimise.";
			label.style.color = "red";
		}
		
		return;
	}
	
	if(sliders.length > datapoints.length)
	{
		var label = document.getElementById("calcstatus");
		if(label)
		{
			label.innerHTML = "The number of data points is smaller than the number of free parameters.";
			label.style.color = "red";
		}
		
		return;
	}
	
	if(datalabels[appmode].length > 0)
	{
		for(co = 0; co < datalabels[appmode].length + (appmode === appmode_receptors ? 1 : 0); co++)
		{
			if(document.getElementById("calcoption" + co).checked) break;
		}
	}
	else
	{
		co = 0;
	}
	
	if(appmode === appmode_receptors && co === 4) co = 5;
	
	var cfg = {
		data: datapoints,
		predict: predict,
		m: m,
		sliderIndices: slider_indices,
		bounds: sliders.map(function(s) {
			return { min: Number(s.min), max: Number(s.max), value: Number(s.value), ext_value: s.ext_value };
		}),
		calcmode: calcmode,
		co: co,
		scale_absolute: scale_absolute,
		appmode: appmode,
		appmode_receptors: appmode_receptors,
		logresid: logresid,
		expparams: expparams,
	};
	
	var fit_result = gridSearch(cfg);
	kp = fit_result.kp;
	resid_min = fit_result.resid_min;
	srsum_min = fit_result.srsum_min;
	var kstep = fit_result.kstep;

	var equal_to_initial = true;
	
	if(calcmode === 3)
	{
		for(i = 0; i < sliders.length; i++)
		{
			if(kp[i] !== Number(sliders[i].value))
			{
				equal_to_initial = false;
				break;
			}
		}
	}
	
	if(calcmode === 4)
	{
		for(i = 0; i < sliders.length; i++)
		{
			if(Math.abs(kp[i] - sliders[i].ext_value) > kstep[i] * 0.1)
			{
				equal_to_initial = false;
				break;
			}
		}
	}
	
	if(calcmode < 4)
	{
		for(i = 0; i < sliders.length; i++)
		{
			sliders[i].value = kp[i].toString();
			slider_input(slider_indices[i], true);
		}
	}
	else
	{
		for(i = 0; i < sliders.length; i++)
		{
			sliders[i].ext_value = kp[i];
		}
	}
	
	var label = document.getElementById("calcstatus");
	
	if(!equal_to_initial) // only in calcmode 3 and 4
	{
		try
		{
			calculate_lsf();
		}
		catch(err)
		{
			// recursion limit reached
			label.innerHTML = "Calculation did not converge or an error occurred.\nTry more optimal starting values.";
			label.style.color = "red";
			
			for(i = 0; i < sliders.length; i++)
			{
				sliders[i].ext_value = undefined;
			}
		}
		
		return;
	}
	
	if(label)
	{
		label.innerHTML = "Calculation finished.";
		label.style.color = "green";
	}
	
	if(calcmode === 4)
	{
		let str = "";
		
		// https://octave.sourceforge.io/optim/function/nlinfit.html
		// https://se.mathworks.com/help/stats/nlinfit.html
		
		// The Jacobian + Gauss-Jordan covariance now lives in core/fit.js. Re-read the
		// (just-written) ext_values into cfg.bounds so the finite-difference deltas match.
		cfg.bounds = sliders.map(function(s) {
			return { min: Number(s.min), max: Number(s.max), value: Number(s.value), ext_value: s.ext_value };
		});
		cfg.srsum_min = srsum_min;
		cfg.resid_min = resid_min;
		var covb = covariance(cfg, kp);
		
		for(i = 0; i < sliders.length; i++)
		{
			if(i > 0) str += "\n\n";
			
			let stdev = covb ? Math.sqrt(covb[i * sliders.length + i]) : NaN;
			let uncer = stdev * tinv_95(resid_min.length - sliders.length);
			let value_str = "";
			let uncer_str = "";
			let rf1, rf2;
			
			if(Number.isFinite(stdev))
			{
				let expstr1 = sliders[i].ext_value.toExponential(0);
				let expstr2 = uncer.toExponential(0);
				
				rf1 = Number(expstr1.substr(expstr1.search(/e/i) + 1));
				rf2 = Number(expstr2.substr(expstr2.search(/e/i) + 1)) - (expstr2.startsWith("1") ? 1 : 0);
				
				while(true)
				{
					value_str = toFixed2(sliders[i].ext_value / Math.pow(10, rf1), rf1 - rf2);
					if(value_str.startsWith("0")) rf1--; else break;
				}
				
				uncer_str = toFixed2(uncer / Math.pow(10, rf1), rf1 - rf2);
			}
			
			str += "Parameter " + document.getElementById("fixlabel" + slider_indices[i]).firstChild.innerHTML +
				"\n\t\u2192 value:  " + format_exponential(sliders[i].ext_value.toExponential(3)) + 
				(Number.isFinite(stdev) ?
					"\n\t\u2192 standard error:  " + format_exponential(stdev.toExponential(3)) + " (" + (100 * stdev / sliders[i].ext_value).toFixed(2) + "%)" + 
					"\n\t\u2192 uncertainty (95%, n=" + (resid_min.length - sliders.length) + "):  " + format_exponential(uncer.toExponential(3)) + " (" + (100 * uncer / sliders[i].ext_value).toFixed(2) + "%)" +
					"\n\t\u2192 rounded:  " + format_exponential("(" + value_str + " \u00B1 " + uncer_str + ")e" + (rf1).toFixed(0))
				: "");
			
			switch(slider_indices[i])
			{
				case 3: S_0 = sliders[i].ext_value; break;
				case 5: E_0 = sliders[i].ext_value; break;
				case 7: K_D = sliders[i].ext_value; break;
				case 9: K_D2 = sliders[i].ext_value; break;
				case 10: Q_0 = sliders[i].ext_value; break;
				default: break;
			}
			
			slider_input(slider_indices[i], true, true);
			
			sliders[i].ext_value = undefined;
		}
		
		alert(str);
	}
}

function alpha(inhib) {
	switch (inhib) {
		case "competitive": return 1 + I0 / Ki;
		case "uncompetitive": return 1;
		case "mixed": return 1 + I0 / Ki;     // simple mixed (α=α′)
		default: return 1;
	}
}

// initial velocity (steady-state MM)
function v0(S) {
	const a = alpha(inhib_mech);
	return (kcat * E0 * S) / (Km * a + S);
}

function mm_rate_curve(x) {      // x = substrate conc
	return v0(x);
}
function ic50_to_Ki(ic50, S, Km) {
	return ic50 / (1 + S / Km);
}
function Ki_to_ic50(Ki, S, Km) {
	return Ki * (1 + S / Km);
}
