// vi: set expandtab shiftwidth=4 softtabstop=4:

/*
 * === UCSF ChimeraX Copyright ===
 * Copyright 2022 Regents of the University of California. All rights reserved.
 * The ChimeraX application is provided pursuant to the ChimeraX license
 * agreement, which covers academic and commercial uses. For more details, see
 * <https://www.rbvi.ucsf.edu/chimerax/docs/licensing.html>
 *
 * This particular file is part of the ChimeraX library. You can also
 * redistribute and/or modify it under the terms of the GNU Lesser General
 * Public License version 2.1 as published by the Free Software Foundation.
 * For more details, see
 * <https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html>
 *
 * THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
 * EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. ADDITIONAL LIABILITY
 * LIMITATIONS ARE DESCRIBED IN THE GNU LESSER GENERAL PUBLIC LICENSE
 * VERSION 2.1
 *
 * This notice must be embedded in or attached to all copies, including partial
 * copies, of the software or any revisions or derivations thereof.
 * === UCSF ChimeraX Copyright ===
 */

// ----------------------------------------------------------------------------
// Routines to interpolate volume data using trilinear interpolation,
// and to interpolate a colormap.  These are for coloring surfaces using
// volume data values.
//
#ifndef INTERPOLATE_HEADER_INCLUDED
#define INTERPOLATE_HEADER_INCLUDED

#include <vector>			// use std::vector<>
#include <arrays/rcarray.h>		// use Numeric_Array

namespace Interpolate
{
enum  Interpolation_Method {INTERP_LINEAR, INTERP_NEAREST};

void interpolate_volume_data(float vertices[][3], int64_t n,
			     float vtransform[3][4],
			     const Reference_Counted_Array::Numeric_Array &data,
			     Interpolation_Method method,
			     float *values, std::vector<int> &outside);

void interpolate_volume_gradient(float vertices[][3], int64_t n,
				 float vtransform[3][4],
				 const Reference_Counted_Array::Numeric_Array &data,
				 Interpolation_Method method,
				 float gradients[][3],
				 std::vector<int> &outside);

void interpolate_colormap(float values[], int64_t n,
			  float color_data_values[], int m,
			  float rgba_colors[][4],
			  float rgba_above_value_range[4],
			  float rgba_below_value_range[4],
			  float rgba[][4]);

void set_outside_volume_colors(int *outside, int64_t n,
			       float rgba_outside_volume[4],
			       float rgba[][4]);

}  // end of namespace interpolate

#endif
