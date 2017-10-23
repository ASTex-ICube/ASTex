
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
