/*******************************************************************************
* ASTex:                                                                       *
* Copyright (C) IGG Group, ICube, University of Strasbourg, France             *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: https://astex-icube.github.io                                      *
* Contact information: astex@icube.unistra.fr                                  *
*                                                                              *
*******************************************************************************/




#include <fstream>
#include <cmath>
#include <algorithm>

#include "distances_maps.h"

namespace ASTex
{

void ASTEX_API normalize_distance_maps(const std::vector< ImageGrayd> input, std::vector<ImageGrayd>& output)
{
	double sum;

	int im_w = input[0].width();
	int im_h = input[0].height();

	int nb_im = input.size();

	output.clear();
	output.resize(nb_im);

	for(int i=0; i< nb_im; ++i)
		output[i].initItk(im_w,im_h,true);

	for (int j = 0; j < im_h; ++j)
	{
		for (int i = 0; i < im_w; ++i)
		{
			sum=0.0;
			for (int k = 0; k < nb_im; ++k)
				sum+= input[k].pixelAbsolute(i,j);

			for (int k = 0; k < nb_im; ++k)
				output[k].pixelAbsolute(i,j)= input[k].pixelAbsolute(i,j)/sum;

		}
	}
}



} //namespace ASTex
