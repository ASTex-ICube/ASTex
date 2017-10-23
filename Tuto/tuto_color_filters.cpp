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
#include <ASTex/colorspace_filters.h>

#include "itkCastImageFilter.h"

using namespace ASTex;


///
/// \brief Example of ColorSpace filters usage
///
/// ColorSpace::FilterXXXtoYYY<TYPE_IMAGE_IN, TYPE_IMAGE_OUT>
///
int main()
{
	ImageRGBd image;
	bool ok = image.load(TEMPO_PATH+"simpleRGB.png");
	if (!ok)
		return 1;

	typedef ImageRGBd::ItkImg IMG_DBL;
	typedef ImageRGBf::ItkImg IMG_FLT;
	typedef ImageRGBu8::ItkImg IMG_U8;

	// first transform [0,255] double -> [0,1] double
	ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::Pointer filter0 =
			ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::New();
	filter0->SetInput(image.itk());


	// RGB double -> XYZ float
	ColorSpace::FilterRGBtoXYZ<IMG_DBL,IMG_FLT>::Pointer filter1 =
			ColorSpace::FilterRGBtoXYZ<IMG_DBL,IMG_FLT>::New();
	filter1->SetInput(filter0->GetOutput());

	// XYZ float -> LUV float
	ColorSpace::FilterXYZtoLUV<IMG_FLT,IMG_FLT>::Pointer filter2 =
			ColorSpace::FilterXYZtoLUV<IMG_FLT,IMG_FLT>::New();
	filter2->SetInput(filter1->GetOutput());

	// LUV float -> XYZ double
	ColorSpace::FilterLUVtoXYZ<IMG_FLT,IMG_DBL>::Pointer filter3 =
			ColorSpace::FilterLUVtoXYZ<IMG_FLT,IMG_DBL>::New();
	filter3->SetInput(filter2->GetOutput());

	// XYZ double -> RGB float
	ColorSpace::FilterXYZtoRGB<IMG_DBL,IMG_FLT>::Pointer filter4 =
			ColorSpace::FilterXYZtoRGB<IMG_DBL,IMG_FLT>::New();
	filter4->SetInput(filter3->GetOutput());

	// finally transform [0,1] float -> [0,255] uint8_t for save
	ColorSpace::FilterRGB01To255<IMG_FLT,IMG_U8>::Pointer filter5 =
			ColorSpace::FilterRGB01To255<IMG_FLT,IMG_U8>::New();
	filter5->SetInput(filter4->GetOutput());


	// create out put ASTex image
	ImageRGBu8 img2(filter5->GetOutput());


	img2.save(TEMPO_PATH+"outFilter.png");

  return EXIT_SUCCESS;
}


