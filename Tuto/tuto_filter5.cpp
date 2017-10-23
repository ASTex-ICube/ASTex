#include <iostream>
#include "itkInPlaceImageFilter.h"
#include "itkSimpleFastMutexLock.h"


#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>
#include <ASTex/exr_io.h>

using namespace ASTex;

/**
 * Example of filter with N (4) outputs
 *  - the 3 first are Grayu8 as second template param of base class (type of first out)
 *  - the 4th is Grayf
 */
class MultipleOutputsFilter: public itk::ImageToImageFilter< ImageRGBu8::ItkImg, ImageGrayu8::ItkImg >
{
public:
	// Standard class typedefs & macros for ikt
	typedef MultipleOutputsFilter                                              Self;
	typedef itk::ImageToImageFilter< ImageRGBu8::ItkImg, ImageGrayu8::ItkImg > Superclass;
	typedef itk::SmartPointer< Self >                                          Pointer;
	typedef itk::SmartPointer< const Self >                                    ConstPointer;
	itkNewMacro(Self);
	itkTypeMacro(MultipleOutputsFilter, ImageToImageFilter);


	/**
	 * @brief convenient function to get output ptr with right type (diff of first out)
	 * @return img ptr with right type
	 */
	ImageGrayf::ItkImg* GetOutput3Grayf()
	{
		return reinterpret_cast<ImageGrayf::ItkImg*>(this->ProcessObject::GetOutput(3));
	}


protected:

	MultipleOutputsFilter()
	{
		this->SetNumberOfRequiredInputs(1);
		this->SetNumberOfRequiredOutputs(4);

		this->SetNthOutput(0,(ImageGrayu8::ItkImg::New()).GetPointer());
		this->SetNthOutput(1,(ImageGrayu8::ItkImg::New()).GetPointer());
		this->SetNthOutput(2,(ImageGrayu8::ItkImg::New()).GetPointer());
		this->SetNthOutput(3,(ImageGrayf::ItkImg::New()).GetPointer());
	}

	virtual ~MultipleOutputsFilter() {}


	//
	// overriden method that generate the output
	//
	void GenerateData() ITK_OVERRIDE
	{
		this->AllocateOutputs();

		ConstImageRGBu8 img_rgb(this->GetInput());
		ImageGrayu8 img_R(this->GetOutput(0));
		ImageGrayu8 img_G(this->GetOutput(1));
		ImageGrayu8 img_B(this->GetOutput(2));
		ImageGrayf img_X(GetOutput3Grayf());

		// apply algo on image
		img_rgb.for_all_pixels([&] (const ImageRGBu8::PixelType& P, int x, int y)
		{
			img_R.pixelAbsolute(x,y) = P[0];
			img_G.pixelAbsolute(x,y) = P[1];
			img_B.pixelAbsolute(x,y) = P[2];
			img_X.pixelAbsolute(x,y) = (float(P[0])+float(P[1])+float(P[2]))/765.0f;
		});
	}

private:
	// to avoid filter object copy
	MultipleOutputsFilter(const Self &);
	void operator=(const Self &);
};



int main()
{

	ImageRGBu8 image;
	bool ok = image.load(TEMPO_PATH+"simpleRGB2.png");
	if (!ok)
		return 1;

	MultipleOutputsFilter::Pointer filter = MultipleOutputsFilter::New();
	// set input of filter
	filter->SetInput(image.itk());

	// create image from output of filter
	ImageGrayu8 img_outR(filter->GetOutput(0));
	ImageGrayu8 img_outG(filter->GetOutput(1));
	ImageGrayu8 img_outB(filter->GetOutput(2));
	ImageGrayf img_outX(filter->GetOutput3Grayf());

	// launch filter
	filter->Update();

	img_outR.save(TEMPO_PATH+"tuto_filter5_R.png");
	img_outG.save(TEMPO_PATH+"tuto_filter5_G.png");
	img_outB.save(TEMPO_PATH+"tuto_filter5_B.png");
	IO::EXR::save(img_outX, TEMPO_PATH+"tuto_filter5_X.exr");

	return EXIT_SUCCESS;
}


