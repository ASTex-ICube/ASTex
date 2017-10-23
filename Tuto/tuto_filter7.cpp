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
#include <ASTex/easy_io.h>

using namespace ASTex;



/**
 * Example of filter with:
 *  - two inputs (A,B) of different types and sizes
 *  - two outputs (C,D) of different types and sizes
 * Template instanciation of ImageToImageFilter is done with types of A & C
 */
class MyFilter: public itk::ImageToImageFilter< ImageRGBu8::ItkImg, ImageRGBd::ItkImg >
{
public:
	// Standard class typedefs & macros for ikt
	typedef MyFilter                                                         Self;
	typedef itk::ImageToImageFilter< ImageRGBu8::ItkImg, ImageRGBd::ItkImg > Superclass;
	typedef itk::SmartPointer< Self >                                        Pointer;
	typedef itk::SmartPointer< const Self >                                  ConstPointer;
	itkNewMacro(Self);
	itkTypeMacro(MyFilter, ImageToImageFilter);


	/// types shortcuts
	using ImA = ConstImageRGBu8::ItkImg;
	using ImB = ConstImageGrayu8::ItkImg;
	using ImC = ImageRGBd::ItkImg;
	using ImD = ImageGrayf::ItkImg;


	void SetInputA(ImA* image)
	{
	 this->SetNthInput(0, image);
	}

	void SetInputB(ImB* image)
	{
	 this->SetNthInput(1, image);
	}

	const ImA* GetInputA()
	{
//		const ImA* ptr = reinterpret_cast<const ImA*>(this->ProcessObject::GetInput(0));
//		return const_cast<ImA*>(ptr);
		return  reinterpret_cast<const ImA*>(this->ProcessObject::GetInput(0));
	}

	const ImB* GetInputB()
	{
//		const ImB* ptr = reinterpret_cast<const ImB*>(this->ProcessObject::GetInput(1));
//		return const_cast<ImB*>(ptr);
		return reinterpret_cast<const ImB*>(this->ProcessObject::GetInput(1));
	}

	ImC* GetOutputC()
	{
		return reinterpret_cast<ImC*>(this->ProcessObject::GetOutput(0));
	}

	ImD* GetOutputD()
	{
		return reinterpret_cast<ImD*>(this->ProcessObject::GetOutput(1));
	}

	void set_size_of_C(int s)
	{
		sz_C_ = s;
	}

	void set_size_of_D(int s)
	{
		sz_D_ = s;
	}


protected:

	itk::SimpleFastMutexLock mutex;

	int sz_C_;
	int sz_D_;

	MyFilter()
	{
		this->SetNumberOfRequiredInputs(2);
		this->SetNumberOfRequiredOutputs(2);
		this->SetNthOutput(0,(ImageRGBd::ItkImg::New()).GetPointer());
		this->SetNthOutput(1,(ImageGrayf::ItkImg::New()).GetPointer());
	}

	virtual ~MyFilter() {}

	// overriden method that generate the output
	void GenerateData() ITK_OVERRIDE
	{
		// get itk pointers
		ImC* itk_outC = GetOutputC();
		ImD* itk_outD = GetOutputD();

		// allocation of output:
		// here not possible to call AllocateOutputs(), we need to allocate ourself
		// C
		Region region;
		region.SetIndex(0,0);
		region.SetIndex(1,0);
		region.SetSize(0,sz_C_);
		region.SetSize(1,sz_C_);
		itk_outC->SetRegions(region);
		itk_outC->Allocate(false);
		// D
		region.SetSize(0,sz_D_);
		region.SetSize(1,sz_D_);
		itk_outD->SetRegions(region);
		itk_outD->Allocate(false);

		// let ITK launch the threads based on the first output cutted in regions.
		// It will compute the first output calling n times ThreadedGenerateData
		this->Superclass::GenerateData();

		// and compute second output after (because size is different)
		ConstImageRGBu8 imgA(GetInputA());
		ConstImageGrayu8 imgB(GetInputB());

		ImageGrayf imgD(GetOutputD());
		imgD.for_all_pixels([&] (float& P, int x, int y)
		{
			const ImageRGBu8::PixelType& Q = imgA.pixelAbsolute(x,y);
			float v = float(imgB.pixelAbsolute(x,y))/255.0f;
			P = (float(Q[0])+float(Q[1])+float(Q[2]))/3.0*v;
		});
	}

	// overriden method that generate the output by regions (of the first output)
	void ThreadedGenerateData(const Region& region, itk::ThreadIdType /*threadId*/) ITK_OVERRIDE
	{
		// create ASTex images
		ConstImageRGBu8 imgA(GetInputA());
		ConstImageGrayu8 imgB(GetInputB());
		ImageRGBd imgC(GetOutputC());

		imgC.for_region_pixels(region, [&] (ImageRGBd::PixelType& P, int x, int y)
		{
			const ImageRGBu8::PixelType& Q = imgA.pixelAbsolute(x,y);
			float v = float(imgB.pixelAbsolute(x,y))/255.0f;
			P[0] = double(Q[0])*v;
			P[1] = double(Q[1])*v;
			P[2] = double(Q[2])*v;
		});
	}

private:
	// to avoid filter object copy
	MyFilter(const Self &);
	void operator=(const Self &);
};



int main()
{

	ImageRGBu8 image;
	bool ok = image.load(TEMPO_PATH+"simpleRGB2.png");
	if (!ok)
		return 1;

	ImageGrayu8 image2(image.width()*2, image.height()*2);
	image2.for_all_pixels([&] (ImageGrayu8::PixelType& Q, int x, int y)
	{
		if (x<y)
			Q = (x%16);
		else
			Q = (y%16);
	});

	MyFilter::Pointer filter = MyFilter::New();
	// set input of filter
	filter->SetInputA(image.itk());
	filter->SetInputB(image2.itk());
	// set param of filter
	filter->set_size_of_C(48);
	filter->set_size_of_D(32);

	// create image from output of filter
	ImageRGBd img_C(filter->GetOutputC());
	ImageGrayf img_D(filter->GetOutputD());
	// launch filter
	filter->Update();

	IO::save01_in_u8(img_C, TEMPO_PATH+"tuto_filter7a.png");
	IO::save01_in_u8(img_D, TEMPO_PATH+"tuto_filter7b.png");

	return EXIT_SUCCESS;
}


