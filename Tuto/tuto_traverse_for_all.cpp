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



#include <iostream>
#include <chrono>

#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>


using namespace ASTex;

// small example to show effect of const usage
double search_max(const ImageGrayd& img)
{
	double maxi = img.pixelAbsolute(0,0);
	// here we have a const ref on image , so parameter of lambda is const ref pixen (or pixel by copy)
	img.for_all_pixels([&maxi] (/*double p*/ const double& p)
	{
		if (p>maxi)
			maxi = p;
	});

	// but this does not compile
	// img.for_all_pixels([&maximums] (double& p,) ...

	return maxi;
}



int main()
{
	ImageGrayd img(64,64);

	// for_all_pixels & for_region_pixels can take:
	// 1 param: the lambda that apply to pixel
	// 2 param: * the lambda that apply to pixel
	//          * a mask (something like (f(x,y)-> bool)

	// the lambda or (functor) can take:
	// 1 param: Pixel& (or const Pixel&)
	// 3 param: Pixel& + x + y (coord. of pixel)


	// traverse all pixels of image with access to each pixel value
	img.for_all_pixels([] (double& p)
	{
		p = 0.0;
	});

	// traverse all pixels of image with access to each pixel value + its coordinates
	img.for_all_pixels([&img] (double& p, int x, int y)
	{
		p = img.width()*y+x;
	});

	// traverse all pixels of image which belong to a mask ( here a lambda)
	img.for_all_pixels([] (double& p)
		{
			p += 0.1;
		},
		[] (int x,int y) { return (x>=0) && (y>=0);} // the mask
	);



	// traverse all pixels of a region of an image
	Region reg = gen_region(4,4,8,8);
	img.for_region_pixels(reg, [] (double& p)
	{
		p += 0.01;
	});

	// or directly
	img.for_region_pixels(4,4,8,8, [] (double& p)
	{
		p += 0.02;
	});


	// parallel traversals:
	//  - whole image
	//  - whole image with mask

	// in parallel mode the lambda or (functor) can take:
	// 1 param: Pixel& (or const Pixel&)
	// 2 param: Pixel& + thread id (0..nbthreads-1)
	// 3 param: Pixel& + x + y
	// 4 param: Pixel& + x + y + thread id

	// parallel traversal
	img.parallel_for_all_pixels([] (double& p)
	{
		p -= 0.1;
	});

	// parallel traversal, with mask (here a lambda)
	img.parallel_for_all_pixels([] (double& p)
	{
		p = 0.0;
	},
	[] (int x,int y) { return (x>=0) && (y>=0);}
	);


	// example: fill with values
	img.parallel_for_all_pixels([] (double& p,int x, int y, uint16_t th)
	{
		p = 2.0 + th + 0.01*y + 0.00001*x;
	});

	// search the max of pixels
	std::vector<double> maximums(std::thread::hardware_concurrency(), std::numeric_limits<double>::min());
	// do not forget to capture maximums by ref in lambda
	img.parallel_for_all_pixels([&maximums] (const double& p, uint16_t th)
	{
		if (p>maximums[th])
			maximums[th] = p; // no pb of concurrency thanks to vector
	});
	for( const auto& m: maximums)
		std::cout << m << "  ";
	// easy way to find max of maximums
	std::sort(maximums.begin(), maximums.end());
	std::cout<< " => " << maximums.back() << std::endl;
	std::cout<< " => " << search_max(img) << std::endl;



	ImageRGBd img2(64,64);

	// extern functions	for_indices allow traversal a region
	// lamdba take 2 param x,y
	// just nice syntax to hide double loop

	// traverse indice [x_min , x_end[ & [y_min , y_end[
	for_indices(0,img.width(),0,img.height(),[&](int i, int j)
	{
		double g = img.pixelAbsolute(i,j);
		img2.pixelAbsolute(i,j) = ImageRGBd::ASTexPixelType(g,g,g);
	});

	// traverse indices of region
	for_region(reg,[&](int i, int j)
	{
		double g = 4.0*img.pixelAbsolute(i,j);
		img2.pixelAbsolute(i,j) = ImageRGBd::ASTexPixelType(g,g,g);
	});




	return EXIT_SUCCESS;
}

