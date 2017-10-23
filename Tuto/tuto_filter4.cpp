#include <iostream>
#include "itkInPlaceImageFilter.h"
#include "itkSimpleFastMutexLock.h"


#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>


using namespace ASTex;


/**
 * Example of filter with two inputs of different types
 * First template param of base class is type of first input
 *
 */
class Diff2InputsFilter: public itk::ImageToImageFilter< ImageRGBu8::ItkImg, ImageGrayu8::ItkImg >
{
public:
	// Standard class typedefs & macros for ikt
	typedef Diff2InputsFilter                                                  Self;
	typedef itk::ImageToImageFilter< ImageRGBu8::ItkImg, ImageGrayu8::ItkImg > Superclass;
	typedef itk::SmartPointer< Self >                                          Pointer;
	typedef itk::SmartPointer< const Self >                                    ConstPointer;
	itkNewMacro(Self);
	itkTypeMacro(Diff2InputsFilter, ImageToImageFilter);


	/**
	 * @brief Set input image A (RGB)
	 * @param image
	 */
	void SetInputImageA(ImageRGBu8::ItkImg* image)
	{
	  this->SetNthInput(0, image);
	}

	/**
	 * @brief Set input image B (Grayd)
	 * @param image
	 */
	void SetInputImageB(ImageGrayd::ItkImg* image)
	{
	  this->SetNthInput(1, image);
	}

	/**
	 * @brief convenient fonction to get pointer on input B (because type is different of input A)
	 * @return pointer of right type
	 */
	const ImageGrayd::ItkImg* GetInputImageB()
	{
		return reinterpret_cast<const ImageGrayd::ItkImg*>(this->ProcessObject::GetInput(1));
	}

protected:

	Diff2InputsFilter()
	{
		this->SetNumberOfRequiredInputs(2);
	}

	virtual ~Diff2InputsFilter() {}


	//
	// overriden method that generate the output
	//
	void GenerateData() ITK_OVERRIDE
	{
		this->AllocateOutputs();

		ConstImageRGBu8 img_A(this->GetInput(0));
		ConstImageGrayd img_B(GetInputImageB());
		ImageGrayu8 img_out(this->GetOutput());

		// apply algo on image
		img_out.for_all_pixels([&] (ImageGrayu8::PixelType& Q, int x, int y)
		{
			ImageRGBu8::PixelType p = img_A.pixelAbsolute(x,y);
			ImageGrayd::PixelType v = img_B.pixelAbsolute(x,y);
			double d = p[0];
			d += p[1];
			d += p[2];
			d/= 3.0;
			Q = d*v;
		});
	}

private:
	// to avoid filter object copy
	Diff2InputsFilter(const Self &);
	void operator=(const Self &);
};



int main()
{

	ImageRGBu8 image;
	bool ok = image.load(TEMPO_PATH+"simpleRGB.png");
	if (!ok)
		return 1;

	ImageGrayd image2(image.width(), image.height());
	int ww = image2.width()/2;
	int hh = image2.height()/2;
	image2.for_all_pixels([&] (ImageGrayd::PixelType& Q, int x, int y)
	{
		double v = (x-ww)*(x-ww) + (y-hh)*(y-hh);
		Q = std::abs(1.0 - std::sqrt(v)/std::max(ww,hh));
	});


	Diff2InputsFilter::Pointer filter = Diff2InputsFilter::New();
	// set input of filter
	filter->SetInputImageA(image.itk());
	filter->SetInputImageB(image2.itk());

	// create image from output of filter
	ImageGrayu8 img_out(filter->GetOutput());
	// launch filter
	filter->Update();

	img_out.save(TEMPO_PATH+"tuto_filter4.png");

	return EXIT_SUCCESS;
}


