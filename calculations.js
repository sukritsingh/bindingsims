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

var datalabels = [];

datalabels[appmode_ligand] = ["[PL]", "[P]"];
datalabels[appmode_homodimer] = ["[P<sub>2</sub>]", "[P]"];
datalabels[appmode_ligands] = ["[PL]", "[PL\u2032]", "[P]"];
datalabels[appmode_receptors] = ["[PL]", "[P\u2032L]", "[P]", "[P\u2032]"];

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
	
	var theor = [];
	
	var slider_indices = [];
	var sliders = [];
	var kmin = [];
	var kmax = [];
	var kstep = [];
	var kp = [];
	var k = [];
	var m = [];
	var fun;
	var fitdata = {};
	
	function fun_ligand_free()
	{
		return calculate_ligand_free(m[5], datapoints[i].x, m[7]);
	}
	
	function fun_ligand_total()
	{
		return calculate_ligand_total(m[5], datapoints[i].x, m[7]);
	}
	
	function fun_homodimer_free()
	{
		return calculate_homodimer_free(datapoints[i].x, m[7]);
	}
	
	function fun_homodimer_total()
	{
		return calculate_homodimer(datapoints[i].x, m[7]);
	}
	
	function fun_ligands()
	{
		return calculate_ligands(datapoints[i].x, m[10], m[3], m[7], m[9]);
	}
	
	function fun_ligands_alt()
	{
		return calculate_ligands(m[5], datapoints[i].x, m[3], m[7], m[9]);
	}
	
	function fun_ligands_free()
	{
		return calculate_ligands_free(datapoints[i].x, m[10], m[3], m[7], m[9]);
	}
	
	switch(appmode)
	{
		case appmode_ligand:
		{
			if(xscale_alternative)
				fun = fun_ligand_free;
			else
				fun = fun_ligand_total;
			break;
		}
		case appmode_homodimer:
		{
			if(xscale_alternative)
				fun = fun_homodimer_free;
			else
				fun = fun_homodimer_total;
			break;
		}
		case appmode_ligands:
		{
			if(xscale_alternative)
				fun = fun_ligands_alt;
			else
				fun = fun_ligands;
			break;
		}
		case appmode_receptors:
		{
			if(xscale_alternative)
				fun = fun_ligands_free;
			else
				fun = fun_ligands;
			break;
		}
		default:
		{
			break;
		}
	}
	
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
	
	for(j = 0; j < (calcmode === 1 ? 2 : 1); j++)
	{
		if(j === 0)
		{
			for(i = 0; i < sliders.length; i++)
			{
				switch(calcmode)
				{
					case 1:
						kmin[i] = Number(sliders[i].min);
						kmax[i] = Number(sliders[i].max);
						kstep[i] = 10;
						break;
					case 2:
						kmin[i] = Number(sliders[i].min);
						kmax[i] = Number(sliders[i].max);
						kstep[i] = 1;
						break;
					case 3:
						kmin[i] = Math.max(Number(sliders[i].min), Number(sliders[i].value) - 5);
						kmax[i] = Math.min(Number(sliders[i].max), Number(sliders[i].value) + 5);
						kstep[i] = 1;
						break;
						
					// extmode only
					case 4:
						kstep[i] = sliders[i].ext_value * 0.0001;
						kmin[i] = sliders[i].ext_value - kstep[i] * 10;
						kmax[i] = sliders[i].ext_value + kstep[i] * 10.1;
						break;
				}
			}
		}
		else
		{
			// calcmode is always 1 in this case
			
			for(i = 0; i < sliders.length; i++)
			{
				kmin[i] = Math.max(Number(sliders[i].min), kp[i] - 20);
				kmax[i] = Math.min(Number(sliders[i].max), kp[i] + 20);
				kstep[i] = 1;
			}
		}
		
		for(i = 0; i < sliders.length; i++)
		{
			k[i] = kmin[i];
		}
		
		var loop = true;
		
		srsum_min = Number.MAX_VALUE;
		resid_min = undefined;
		
		while(loop)
		{
			for(i = 0; i < sliders.length; i++)
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
			
			for(i = 0; i < datapoints.length; i++)
			{
				var cd = fun();
				var divisor;
				
				if(scale_absolute || appmode === appmode_receptors)
					divisor = 1;
				else
					divisor = cd.d.reduce(function(t,v,i){ return t + v * cd.m[i]; }, 0);
				
				theor[i] = {
					x: datapoints[i].x,
					y: cd.d[co] / divisor * (scale_absolute || appmode === appmode_receptors ? 1 : cd.m[co])
				};
			}
			
			fitdata = sum_squared_residuals(datapoints, theor, logresid);
			
			if(fitdata.sum < srsum_min)
			{
				srsum_min = fitdata.sum;
				resid_min = fitdata.residuals;
				kp = k.slice();
			}
			
			for(i = 0; i < sliders.length; i++)
			{
				if((k[i] += kstep[i]) <= kmax[i])
					break;
				else if(i < sliders.length - 1)
					k[i] = kmin[i];
				else
					loop = false;
			}
		}
	}
	
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
		
		var covb;
		
		if(resid_min.length - sliders.length > 0)
		{
			// calculate mean square error (mse)
			
			let mse = srsum_min / (resid_min.length - sliders.length);
			
			// calculate jacobian matrix J
			
			let J_rows = resid_min.length;
			let J_cols = sliders.length;
			let J = new Array(J_rows * J_cols);
			
			for(i = 0; i < sliders.length; i++)
			{
				m[slider_indices[i]] = kp[i];
			}
			
			for(i = 0; i < datapoints.length; i++)
			{
				for(j = 0; j < sliders.length; j++)
				{
					let delta = sliders[j].ext_value * 0.0001;
					let temp = m[slider_indices[j]];
					
					m[slider_indices[j]] = temp - delta;
					var cd1 = fun();
					
					m[slider_indices[j]] = temp + delta;
					var cd2 = fun();
					
					m[slider_indices[j]] = temp;
					
					var divisor1, divisor2;
					
					if(scale_absolute || appmode === appmode_receptors)
					{
						divisor1 = 1;
						divisor2 = 1;
					}
					else
					{
						divisor1 = cd1.d.reduce(function(t,v,i){ return t + v * cd1.m[i]; }, 0);
						divisor2 = cd2.d.reduce(function(t,v,i){ return t + v * cd2.m[i]; }, 0);
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
			
			/*
			console.log(mse);
			console.log(J);
			console.log(JJ);
			console.log(covb);
			*/
		}
		
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
