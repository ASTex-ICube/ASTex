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
#include "itkInPlaceImageFilter.h"
#include "itkSimpleFastMutexLock.h"


#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>


using namespace ASTex;


/**
 * Example of filter with two inputs
 *   RGBu8 + RGBu8 -> Grayu8
 */
class TwoInputsFilter: public itk::ImageToImageFilter< ImageRGBu8::ItkImg, ImageGrayu8::ItkImg >
{
public:
	// Standard class typedefs & macros for ikt
	typedef TwoInputsFilter                                                    Self;
	typedef itk::ImageToImageFilter< ImageRGBu8::ItkImg, ImageGrayu8::ItkImg > Superclass;
	typedef itk::SmartPointer< Self >                                          Pointer;
	typedef itk::SmartPointer< const Self >                                    ConstPointer;
	itkNewMacro(Self);
	itkTypeMacro(TwoInputsFilter, ImageToImageFilter);


	void SetInputA(ImageRGBu8::ItkImg* image)
	{
	  this->SetNthInput(0, image);
	}

	void SetInputB(ImageRGBu8::ItkImg* image)
	{
	  this->SetNthInput(1, image);
	}

protected:

	// protected constructor (forbid usage of new and variable declaration)
	TwoInputsFilter()
	{
		this->SetNumberOfRequiredInputs(2);
	}

	virtual ~TwoInputsFilter() {}


	//
	// overriden method that generate the output
	//
	void GenerateData() ITK_OVERRIDE
	{
		this->AllocateOutputs();

		ConstImageRGBu8 img_A(this->GetInput(0));
		ConstImageRGBu8 img_B(this->GetInput(1));
		ImageGrayu8 img_out(this->GetOutput());

		// apply algo on region
		img_out.for_all_pixels([&] (ImageGrayu8::PixelType& Q, int x, int y)
		{
			auto Pa = img_A.pixelEigenAbsolute(x,y);
			auto Pb = img_B.pixelEigenAbsolute(x,y);
			Pa /= 2;
			Pb /= 2;
			auto p = Pa+Pb;
			double d = (p[0] + p[1] + p[2])/3.0;
			Q = ImageGrayu8::itkPixel(d);
		});
	}

private:
	// to avoid filter object copy
	TwoInputsFilter(const Self &);
	void operator=(const Self &);
};



int main()
{

	ImageRGBu8 image;
	bool ok = image.load(TEMPO_PATH+"simpleRGB.png");
	if (!ok)
		return 1;

	ImageRGBu8 image2;
	ok = image2.load(TEMPO_PATH+"simpleRGB2.png");
	if (!ok)
		return 1;

	TwoInputsFilter::Pointer filter = TwoInputsFilter::New();
	// set input of filter
	filter->SetInputA(image.itk());
	filter->SetInputB(image2.itk());

	// create image from output of filter
	ImageGrayu8 img_out(filter->GetOutput());

	// launch filter
	filter->Update();

	img_out.save(TEMPO_PATH+"tuto_filter3.png");

	return EXIT_SUCCESS;
}


