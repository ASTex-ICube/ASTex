#include <iostream>
#include "itkInPlaceImageFilter.h"
#include "itkSimpleFastMutexLock.h"


#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>
#include <ASTex/exr_io.h>

using namespace ASTex;


/**
 * Example of Filter with an output with size different from input (see allocation in GenerateData)
 */
class DiffSizeOutputFilter: public itk::ImageToImageFilter< ImageRGBu8::ItkImg, ImageGrayu8::ItkImg >
{
public:
	// Standard class typedefs & macros for ikt
	typedef DiffSizeOutputFilter                                                     Self;
	typedef itk::ImageToImageFilter< ImageRGBu8::ItkImg, ImageGrayu8::ItkImg > Superclass;
	typedef itk::SmartPointer< Self >                                          Pointer;
	typedef itk::SmartPointer< const Self >                                    ConstPointer;
	itkNewMacro(Self);
	itkTypeMacro(DiffSizeOutputFilter, ImageToImageFilter);

protected:

	DiffSizeOutputFilter()
	{}

	virtual ~DiffSizeOutputFilter() {}

	void GenerateData() ITK_OVERRIDE
	{
		// get itk pointers
		ConstImageRGBu8 img_in(this->GetInput());
		ImageGrayu8 img_out(this->GetOutput(0));

		// allocation of output:
		// here not possible to call AllocateOutputs(), we need to allocate ourself
		Region region;
		region.SetIndex(0,0);
		region.SetIndex(1,0);
		region.SetSize(0,img_in.width()/2);
		region.SetSize(1,img_in.height()/2);
		img_out.itk()->SetRegions(region);
		img_out.itk()->Allocate(false);

		// apply algo on image
		img_out.for_all_pixels([&] (ImageGrayu8::PixelType& P, int x, int y)
		{
			const ImageRGBu8::PixelType& Q1 = img_in.pixelAbsolute(2*x,2*y);
			const ImageRGBu8::PixelType& Q2 = img_in.pixelAbsolute(2*x+1,2*y);
			const ImageRGBu8::PixelType& Q3 = img_in.pixelAbsolute(2*x,2*y+1);
			const ImageRGBu8::PixelType& Q4 = img_in.pixelAbsolute(2*x+1,2*y+1);
			float P1 = float(Q1[0])+float(Q1[1])+float(Q1[2]);
			float P2 = float(Q2[0])+float(Q2[1])+float(Q2[2]);
			float P3 = float(Q3[0])+float(Q3[1])+float(Q3[2]);
			float P4 = float(Q4[0])+float(Q4[1])+float(Q4[2]);
			P = (P1+P2+P3+P4)/12.0f;
		});
	}

private:
	// to avoid filter object copy
	DiffSizeOutputFilter(const Self &);
	void operator=(const Self &);
};



/**
 * Example of MT-Filter with an output with size different from input (see allocation in GenerateData)
 */
class DiffSizeOutputMTFilter: public itk::ImageToImageFilter< ImageRGBu8::ItkImg, ImageGrayu8::ItkImg >
{
public:
	// Standard class typedefs & macros for ikt
	typedef DiffSizeOutputMTFilter                                                     Self;
	typedef itk::ImageToImageFilter< ImageRGBu8::ItkImg, ImageGrayu8::ItkImg > Superclass;
	typedef itk::SmartPointer< Self >                                          Pointer;
	typedef itk::SmartPointer< const Self >                                    ConstPointer;
	itkNewMacro(Self);
	itkTypeMacro(DiffSizeOutputMTFilter, ImageToImageFilter);

protected:

	DiffSizeOutputMTFilter()
	{}

	virtual ~DiffSizeOutputMTFilter() {}

	// overriden method that generate the output
	void GenerateData() ITK_OVERRIDE
	{
		// get itk pointers
		const ImageRGBu8::ItkImg* itk_in = this->GetInput();
		ImageGrayu8::ItkImg* itk_out =this->GetOutput();

		// allocation of output:
		// here not possible to call AllocateOutputs(), we need to allocate ourself
		Region region;
		region.SetIndex(0,0);
		region.SetIndex(1,0);
		region.SetSize(0,itk_in->GetLargestPossibleRegion().GetSize()[0]/2);
		region.SetSize(1,itk_in->GetLargestPossibleRegion().GetSize()[1]/2);
		itk_out->SetRegions(region);
		itk_out->Allocate(false);

		// let ITK launch the thread (see just bellow)
		Superclass::GenerateData();
	}

	// overriden method that generate the output by region
	void ThreadedGenerateData(const Region& region, itk::ThreadIdType /*threadId*/) ITK_OVERRIDE
	{
		// create ASTex images
		ConstImageRGBu8 img_in(this->GetInput());
		ImageGrayu8 img_out(this->GetOutput(0));

		// apply algo on region
		img_out.for_region_pixels(region, [&] (ImageGrayu8::PixelType& P, int x, int y)
		{
			const ImageRGBu8::PixelType& Q1 = img_in.pixelAbsolute(2*x,2*y);
			const ImageRGBu8::PixelType& Q2 = img_in.pixelAbsolute(2*x+1,2*y);
			const ImageRGBu8::PixelType& Q3 = img_in.pixelAbsolute(2*x,2*y+1);
			const ImageRGBu8::PixelType& Q4 = img_in.pixelAbsolute(2*x+1,2*y+1);
			float P1 = float(Q1[0])+float(Q1[1])+float(Q1[2]);
			float P2 = float(Q2[0])+float(Q2[1])+float(Q2[2]);
			float P3 = float(Q3[0])+float(Q3[1])+float(Q3[2]);
			float P4 = float(Q4[0])+float(Q4[1])+float(Q4[2]);
			P = (P1+P2+P3+P4)/12.0f;
		});
	}

private:
	// to avoid filter object copy
	DiffSizeOutputMTFilter(const Self &);
	void operator=(const Self &);
};



int main()
{

	ImageRGBu8 image;
	bool ok = image.load(TEMPO_PATH+"simpleRGB2.png");
	if (!ok)
		return 1;

	DiffSizeOutputFilter::Pointer filter = DiffSizeOutputFilter::New();
	// set input of filter
	filter->SetInput(image.itk());

	// create image from output of filter
	ImageGrayu8 img_out(filter->GetOutput());

	// launch filter
	filter->Update();

	img_out.save(TEMPO_PATH+"tuto_filter6a.png");


	//with MT Filter
	DiffSizeOutputMTFilter::Pointer filtermt = DiffSizeOutputMTFilter::New();

	// set input of filter
	filtermt->SetInput(image.itk());

	// create image from output of filter
	ImageGrayu8 img_out2(filtermt->GetOutput());

	// launch filter
	filtermt->Update();

	img_out2.save(TEMPO_PATH+"tuto_filter6b.png");



	return EXIT_SUCCESS;
}


