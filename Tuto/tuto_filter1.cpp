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
#include "itkImageToImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithOnlyIndex.h"
#include "itkSimpleFastMutexLock.h"

#include <itkMetaDataDictionary.h>
#include <itkMetaDataObject.h>


#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>


using namespace ASTex;


/**
 * Example of ITK Filter
 *  RGBu8 -> Grayu8
 */
class MyImageFilter: public itk::ImageToImageFilter< ImageRGBu8::ItkImg, ImageGrayu8::ItkImg >
{
public:
	// Filter is not template so define input & output types for more simple writing
	using TInputImage  = ImageRGBu8::ItkImg;
	using TOutputImage = ImageGrayu8::ItkImg;
	using TInputPixel  = TInputImage::PixelType;
	using TOutputPixel = TOutputImage::PixelType;

	// Standard class typedefs & macros for ikt
	// this is necessary for filter usage
	typedef MyImageFilter                                          Self;
	typedef itk::ImageToImageFilter< TInputImage, TOutputImage >   Superclass;
	typedef itk::SmartPointer< Self >                              Pointer;
	typedef itk::SmartPointer< const Self >                        ConstPointer;
	itkNewMacro(Self);
	itkTypeMacro(MyImageFilter, ImageToImageFilter);

	/**
	 * @brief set center for filter computation
	 * @param c
	 */
	void setCenter(const itk::Offset<2>& c)
	{
		center_ = c;
	}

	/**
	 * @brief set radius for filter computation
	 * @param r
	 */
	void setRadius(long r)
	{
		radius_ = r;
	}

protected:
	/// filter data
	itk::Offset<2> center_;
	long radius_;


	/// protected constructor (forbid usage of new and variable declaration)
	MyImageFilter():
		center_({{0,0}}),
		radius_(1)
	{}

	virtual ~MyImageFilter() {}


	//
	// overriden method that generate the output
	//
	void GenerateData() ITK_OVERRIDE
	{
		//  we are responsible of allocating the images buffers
		// in this simple case (output of same size than inputs) just call:
		this->AllocateOutputs();

		// get in & out images ptr
		TInputImage::ConstPointer input  = this->GetInput();
		TOutputImage::Pointer output     = this->GetOutput();

		// define in & out iterator on image
		itk::ImageRegionConstIterator< TInputImage >  itIn  (input, input->GetBufferedRegion());
		itk::ImageRegionIterator< TOutputImage > itOut (output, output->GetBufferedRegion());

		// init
		itIn.GoToBegin();
		itOut.GoToBegin();

		// traverse image
		while (!itOut.IsAtEnd())
		{
			// get input pixel
			const TInputPixel& p = itIn.Get();

			// some simple computation
			double d = p[0];
			d += p[1];
			d += p[2];
			d/= 3.0;

			// get index (position of pixel)
			// warning not efficient here because with RegionIterator it is computed
			// use RegionIteratorWithIndex instead
			itk::Index<2> pos = itIn.GetIndex();

			// fading
			pos -= center_;
			double dec = 1.0 - sqrt(pos[0]*pos[0] + pos[1]*pos[1])/radius_;
			d *= dec;

			// compute gray pixel
			TOutputPixel q(d);

			// set to ouput image
			itOut.Set(q);

			// next pixel
			++itIn;
			++itOut;
		}
	}

private:
	// to avoid filter object copy
	MyImageFilter(const Self &);
	void operator=(const Self &);
};



/**
 * Example of Multi-threaded ITK Filter
 *   RGBu8 -> Grayu8
 */
class MyMTImageFilter: public itk::ImageToImageFilter< ImageRGBu8::ItkImg, ImageGrayu8::ItkImg >
{
public:
	// Filter is not template so define input & output types for more simple writing
	using TInputImage  = ImageRGBu8::ItkImg;
	using TOutputImage = ImageGrayu8::ItkImg;
	using TInputPixel  = TInputImage::PixelType;
	using TOutputPixel = TOutputImage::PixelType;

	// Standard class typedefs & macros for ikt
	typedef MyMTImageFilter                                        Self;
	typedef itk::ImageToImageFilter< TInputImage, TOutputImage >   Superclass;
	typedef itk::SmartPointer< Self >                              Pointer;
	typedef itk::SmartPointer< const Self >                        ConstPointer;
	itkNewMacro(Self);
	itkTypeMacro(MyImageFilter, ImageToImageFilter);

	/**
	 * @brief set center for filter computation
	 * @param c
	 */
	void setCenter(const itk::Offset<2>& c)
	{
		center_ = c;
	}

	/**
	 * @brief set radius for filter computation
	 * @param r
	 */
	void setRadius(long r)
	{
		radius_ = r;
	}

protected:
	/// filter data
	itk::Offset<2> center_;
	long radius_;

	/// just for nice cout
	itk::SimpleFastMutexLock mutex;


	/// protected constructor (forbid usage of new and variable declaration)
	MyMTImageFilter():
		center_({{0,0}}),
		radius_(1)
	{}

	virtual ~MyMTImageFilter() {}


	//
	// overriden method that generate the output (by region)
	// if GenerateData is not overriden it calls AllocateOutputs, then
	// it cuts the (first) output into nbthreads regions and call ThreadedGenerateData
	// for each of them, in // of course
	//
	void ThreadedGenerateData(const Region& region, itk::ThreadIdType threadId) ITK_OVERRIDE
	{
		mutex.Lock();
		std::cout << "Thread " << threadId << " given region: " << region << std::endl;
		mutex.Unlock();

		// get in & out images ptr
		TInputImage::ConstPointer input  = this->GetInput();
		TOutputImage::Pointer output     = this->GetOutput();

		// define in & out iterator on image
		itk::ImageRegionConstIterator< TInputImage >  itIn  (input, region);
		itk::ImageRegionIterator< TOutputImage > itOut (output, region);

		// init
		itIn.GoToBegin();
		itOut.GoToBegin();

		// traverse region
		while (!itOut.IsAtEnd())
		{
			// get input pixel
			const TInputPixel& p = itIn.Get();

			// some simple computation
			double d = p[0];
			d += p[1];
			d += p[2];
			d/= 3.0;

			// get index (position of pixel)
			// warning not efficient here because with RegionIterator it if computed
			// use RegionIteratorWithIndex instead
			itk::Index<2> pos = itIn.GetIndex();

			// fading
			pos -= center_;
			double dec = std::abs(1.0 - std::sqrt(pos[0]*pos[0] + pos[1]*pos[1])/radius_);
			d *= dec;

			// compute gray pixel
			TOutputPixel q(d);

			// set to ouput image
			itOut.Set(q);

			// next pixel
			++itIn;
			++itOut;
		}
	}

private:
	// to avoid filter object copy
	MyMTImageFilter(const Self &);
	void operator=(const Self &);
};




int main()
{

	ImageRGBu8 image;

	bool ok = image.load(TEMPO_PATH+"simpleRGB.png");
	if (!ok)
		return 1;

	long w = image.width();
	long h = image.height();

	MyImageFilter::Pointer filter = MyImageFilter::New();

	// set input of filter
	filter->SetInput(image.itk());
	// add external info for filter
	filter->setCenter(Offset({w/2,h/2}));
	// add external info
	filter->setRadius(w/1.414);

	// create image from output of filter
	ImageGrayu8 img2(filter->GetOutput());
	// launch filter
	filter->Update();

	img2.save(TEMPO_PATH+"tuto_filter1a.png");

	// with MT Filter

	MyMTImageFilter::Pointer filtermt = MyMTImageFilter::New();

	// set input of filter
	filtermt->SetInput(image.itk());
	// add external info for filter
	filtermt->setCenter(Offset({w/3,h/3}));
	// add external info
	filtermt->setRadius(w/2.828);

	// create image from output of filter
	ImageGrayu8 img3(filtermt->GetOutput());
	// launch filter
	filtermt->Update();

	img3.save(TEMPO_PATH+"tuto_filter1b.png");


	return EXIT_SUCCESS;
}


