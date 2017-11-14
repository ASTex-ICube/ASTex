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
#include <ASTex/image_rgb.h>

using namespace ASTex;


/**
 * usage of iterators for const traversals (no modif of image)
 * Define where to iterate:
 *  -  it = image.beginConstIterator()
 *  -  it = image.beginConstIterator(x,y,w,h)
 *  -  it = image.beginConstIterator(reg)
 * Test for ending:
 *  -  !it.IsAtEnd()
 * Iter: ++it
 */
void test_const_iterator(const ImageRGBu8& image)
{
	uint32_t nb = 0;
	auto average = ImageRGBu8::eigenPixel(0);

	// The three following usages of ConstIterator vary only with parameter of beginConstIterator
	// no param, 4 int, or Region

	// 1 - IMAGE
	// traverse whole image with iterator
	// we have to use ConstIterator because image is const &
	for (ImageRGBu8::ConstIterator it = image.beginConstIterator(); !it.IsAtEnd(); ++it)
	{
		// construct a RGBd from itValue for using +=
		average += ImageRGBu8::eigenPixel(it.Value());
		nb ++;
	}
	average /= nb;
	std::cout << "Average: "<<ImageRGBu8::itkPixel(average) << std::endl;

	uint32_t W = image.width();
	uint32_t H = image.height();
	nb = 0;
	average = ImageRGBu8::eigenPixel(0);

	// 2 - ZONE x,y,w,h
	// traverse a part of image (a quad of size half of image centered)
	// Region is given by  begin pos x, begin pos y, width, height.
	for (ImageRGBu8::ConstIterator it = image.beginConstIterator(W/4, H/4, W/2, H/2); !it.IsAtEnd(); ++it)
	{
		average += ImageRGBu8::eigenPixel(it.Value());
		nb ++;
	}
	average /= nb;
	std::cout << "Average; "<< ImageRGBu8::itkPixel(average) << std::endl;


	Region reg= gen_region(W/4, H/4, W/2, H/2);
	nb = 0;
	average = ImageRGBu8::eigenPixel(0);

	// 3 - Region
	// traverse a region
	for (ImageRGBu8::ConstIterator it = image.beginConstIterator(reg); !it.IsAtEnd(); ++it)
	{
		average += ImageRGBu8::eigenPixel(it.Value());
		nb ++;
	}
	average /= nb;
	std::cout << "Average; "<< ImageRGBu8::itkPixel(average) << std::endl;

}

/**
 * usage of iterator for modifications
 */
void test_iterator(ImageRGBu8& image)
{
	// traverse image with modification: useIterator
	for (ImageRGBu8::Iterator it = image.beginIterator(); !it.IsAtEnd(); ++it)
	{
		it.Value() = ImageRGBu8::itkPixel(ImageRGBu8::eigenPixel(it.Value())/2);
	}
	// Same possible parameter than with ConstIterators (x,y,w,h) (Region)
}


/**
 * usage of several iterators concurrently
 */
void test_many_iterators(const ImageRGBu8& input1, const ImageRGBu8& input2, ImageRGBd& output)
{
	// check size compatibility
	if (!(output.is_initialized_as(input1))&&(output.is_initialized_as(input2)))
		return;

	// traverse the three images
	auto it_in1 = input1.beginConstIterator();
	auto it_in2 = input2.beginConstIterator();
	auto it_out = output.beginIterator();

	while (!it_out.IsAtEnd())
	{
		it_out.Value() = ImageRGBu8::itkPixel((ImageRGBu8::eigenPixel(it_in1.Value()) + ImageRGBu8::eigenPixel(it_in2.Value()))/512);
		 ++it_out,++it_in1,++it_in2;
	}

	// or, for those who do not like to pollute global scope, use for & tuple:
	// 0 : output, 1 input1, 2 input2
	for (auto its = std::make_tuple(output.beginIterator(),input1.beginConstIterator(),input2.beginConstIterator());
		 !std::get<0>(its).IsAtEnd(); ++std::get<0>(its), ++std::get<1>(its), ++std::get<2>(its))
	{
		std::get<0>(its).Value() = ImageRGBu8::itkPixel((ImageRGBu8::eigenPixel(std::get<1>(its).Value()) + ImageRGBu8::eigenPixel(std::get<2>(its).Value()))/512);


	}

}


int main()
{
	ImageRGBu8 image;

	bool ok = image.load(TEMPO_PATH+"simpleRGB.png");
	if (!ok)
		return 1;

	test_const_iterator(image);
	test_iterator(image);
	test_const_iterator(image);

	ImageRGBu8 image2;
	image2.share(image); // same data

	ImageRGBd image3;
	test_many_iterators(image,image2, image3);

  return EXIT_SUCCESS;
}



