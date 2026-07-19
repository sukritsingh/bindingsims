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

// core/format.js — number/exponential/superscript formatting and the rich-text
// SVG renderer.
//
// Moved verbatim from simulation.js / calculations.js in refactor/01-core-math;
// the declarations below are byte-identical. `render_text_svg` still reads
// `svg_xmlns`, `delete_all_children` and `document` from the shared global scope
// (declared in simulation.js / the page), exactly as it did before the move.

function render_text_svg(ele, str)
{
	var fontfamily = "";
	var normalfont = "12px";
	var scriptfont = "10px";
	
	var supbl = -4;
	var subpl = 2.5;
	
	var nextchar = 0;
	var prevchar = -1;
	var search = 0;
	var token = "";
	var baseline = 0;
	var c;
	
	delete_all_children(ele);
	
	do
	{
		search = str.substr(prevchar + 1).search(/[\uEEE0-\uEEEF]/);
		nextchar = prevchar + 1 + search;
		token = str.substring(prevchar + 1, (search === -1) ? undefined : nextchar);
		
		c = (prevchar === -1 ? 0xEEE0 : str.charCodeAt(prevchar));
		
		function new_textspan(bold, italic, script, align, content)
		{
			if(!content) return;

			var span = document.createElementNS(svg_xmlns, "tspan");
			if(bold) span.style.fontWeight = "bold";
			if(italic) span.style.fontStyle = "italic";
			if(script) span.style.fontSize = "10px";
			if(align - baseline !== 0) span.setAttribute("dy", align - baseline);
			baseline = align;

			var textnode = document.createTextNode(content);
			span.appendChild(textnode);

			ele.appendChild(span);
		}

		switch(c)
		{
			case 0xEEE0: // normal
			{
				new_textspan(false, false, false, 0, token);
				break;
			}
			case 0xEEE1: // superscript
			{
				new_textspan(false, false, true, supbl, token);
				break;
			}
			case 0xEEE2: // subscript
			{
				new_textspan(false, false, true, subpl, token);
				break;
			}
			case 0xEEE4: // italic
			{
				new_textspan(false, true, false, 0, token);
				break;
			}
			case 0xEEE5: // italic superscript
			{
				new_textspan(false, true, true, supbl, token);
				break;
			}
			case 0xEEE6: // italic subscript
			{
				new_textspan(false, true, true, subpl, token);
				break;
			}
			case 0xEEE8: // bold normal
			{
				new_textspan(true, false, false, 0, token);
				break;
			}
			case 0xEEE9: // bold superscript
			{
				new_textspan(true, false, true, supbl, token);
				break;
			}
			case 0xEEEA: // bold subscript
			{
				new_textspan(true, false, true, subpl, token);
				break;
			}
			case 0xEEEC: // bold italic
			{
				new_textspan(true, true, false, 0, token);
				break;
			}
			case 0xEEED: // bold italic superscript
			{
				new_textspan(true, true, true, supbl, token);
				break;
			}
			case 0xEEEE: // bold italic subscript
			{
				new_textspan(true, true, true, subpl, token);
				break;
			}
			default:
			{
				prevchar = nextchar;
				continue;
			}
		}
		
		prevchar = nextchar;
	}
	while(search !== -1);
}

function replace_minus_signs(str)
{
	return str.replace(/\u002D/g, "\u2212");
}

function number_to_superscript(str)
{
	return str.replace(/\u002D/g, "\u207B").replace(/\u002B/g, "\u207A").replace(/\u0030/g, "\u2070").replace(/\u0031/g, "\u00B9").replace(/\u0032/g, "\u00B2").replace(/\u0033/g, "\u00B3").replace(/\u0034/g, "\u2074").replace(/\u0035/g, "\u2075").replace(/\u0036/g, "\u2076").replace(/\u0037/g, "\u2077").replace(/\u0038/g, "\u2078").replace(/\u0039/g, "\u2079");
}

function toFixed2(num, digits)
{
	if(digits >= 0)
	{
		return num.toFixed(digits);
	}
	else
	{
		let d = Math.pow(10, Math.max(digits, -Math.floor(Math.log10(num))));
		return (Math.round(num * d) / d).toFixed(0);
	}
}

function zeropad(number, digits)
{
	var str = number.toFixed(0);
	return "0".repeat(str.length < digits ? digits - str.length : 0) + str;
}

function format_exponential(str)
{
	let epos = str.search(/e/i);
	let exp = Number(str.substring(epos + 1));
	
	if(exp == 0)
		if(str.charAt(0) === "(" && str.charAt(epos - 1) === ")")
			return str.substring(1, epos - 1);
		else
			return str.substring(0, epos);
	else
		return str.substring(0, epos) + " \u22C5 10" + number_to_superscript((exp).toFixed(0));
}

export {
	render_text_svg, replace_minus_signs, number_to_superscript,
	toFixed2, zeropad, format_exponential,
};
