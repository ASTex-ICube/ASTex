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



int main()
{
	// All following explanation could also apply to RGBA

	using IMG = ImageRGBu8;
	using PIX = IMG::PixelType;
	using DPE = IMG::DoublePixelEigen;

	IMG image(4,4);
	// PIX => itk::RGBPixel<uint8>

	// pixelAbsolute return ref or const-ref to PixelType
	PIX& p1 = image.pixelAbsolute(0,0);
	PIX& p2 = image.pixelAbsolute(1,0);
	PIX& p3 = image.pixelAbsolute(2,0);
	PIX& p4 = image.pixelAbsolute(3,0);

	ImageGrayu8::PixelType pg;


//	auto Q1 =eigenPixel<double>(pg);

	// iterator.Value() return (const) PixelType
	IMG::ConstIterator it = image.beginConstIterator();
	const PIX& v=it.Value();
	std::cout << v << std::endl;

	// for_all traversal fonction lambda take PixelType parameter for reference (const)

	auto pi = IMG::itkPixel(65, 55, 45);
	image.for_all_pixels([&] (PIX& p)
	{
		p = pi;
		pi[0] += 5;
		pi[1] += 15;
		pi[2] += 25;
	});

	// Examples of function-constructors
	p1 = IMG::itkPixel(134);
	p2 = IMG::itkPixel(128,128,128);

	// Why using DoublePixelEigen / LongPixelEigen type

	p3 = (p1 + p2); // possible overflow !
//	p3 /= 2; not possible with itk pixel type
	std::cout<< p1 << " + "<< p2 << " (compute uint8) = " << p3 << std::endl;

	// easy type conversion:
	p3 = IMG::itkPixel((eigenPixel<int64_t>(p1)+eigenPixel<int64_t>(p2))/2);
	p4 = IMG::itkPixel((eigenPixel<double>(p1)+eigenPixel<double>(p2))/2);

	std::cout<< p1 << " + "<< p2 << " /2 (compute eigen) = " << p3 << " / " << p4 << std::endl;

	// easy direct access with on the fly conversion
	DPE dp1 = image.pixelEigenAbsolute(0,0);
	DPE dp2 = image.pixelEigenAbsolute(1,0); // must be cast or affected before used (proxy)
	DPE dp3 = image.pixelEigenAbsolute(2,0);
	image.pixelEigenAbsolute(1,1) = (dp1+dp2+dp3)/3;
	std::cout<< image.pixelAbsolute(0,0) << " + "<< image.pixelAbsolute(1,0) << " + "<< image.pixelAbsolute(2,0) << " /3 = " << image.pixelAbsolute(1,1) << std::endl;


	auto readEigen = [&image] (int i,int j)	{return DPE(image.pixelEigenAbsolute(i,j));};
	auto read01Eigen = [&image] (int i,int j){return DPE(image.pixelNormEigenAbsolute(i,j));};

	// easy on the fly conversion to 0-1 (use transpose for horizontal output)
	std::cout << read01Eigen(0, 0).transpose() << " + " << read01Eigen(1, 0).transpose() << " + " << read01Eigen(2, 0).transpose() << " /3 = " << read01Eigen(1, 1).transpose() << std::endl;

	image.pixelNormEigenAbsolute(2, 2) = DPE(read01Eigen(0, 0)/2.0);
	std::cout << read01Eigen(2, 2).transpose() << " -> " << image.pixelAbsolute(2, 2) << std::endl;

	// work with all type of image
	ImageRGBd imaged(4, 4);
	imaged.pixelEigenAbsolute(1, 1) = (readEigen(0, 0) + readEigen(1, 0) + readEigen(2, 0)) / 3;

	return EXIT_SUCCESS;
}

