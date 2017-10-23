#include <cmath>

#include "fourier.h"
#include "region_traversor.h"

// IMAGE
#include <itkShiftScaleImageFilter.h>
#include <itkAddImageFilter.h>
//#include <itkMultiplyImageFilter.h>
#include <itkSquareImageFilter.h>
#include <itkSqrtImageFilter.h>
#include <itkRegionOfInterestImageFilter.h>
#include <itkStatisticsImageFilter.h>

// FFT
#include <itkForwardFFTImageFilter.h>
#include <itkInverseFFTImageFilter.h>
#include <itkFFTShiftImageFilter.h>

// COMPLEX
#include <itkComplexToRealImageFilter.h>
#include <itkComplexToImaginaryImageFilter.h>
#include <itkComplexToModulusImageFilter.h>
#include <itkComplexToPhaseImageFilter.h>
#include <itkMagnitudeAndPhaseToComplexImageFilter.h>
#include <itkForwardFFTImageFilter.h>
#include <itkInverseFFTImageFilter.h>
#include <itkFFTShiftImageFilter.h>

namespace ASTex
{

namespace Fourier
{

void ASTEX_API autoCorrelation_normalize(ImageGrayd& acorr, double v)
{
	assert (v!=0);
	acorr.for_all_pixels([&](double& P)
	{
		P*=v;
	});
}

void ASTEX_API autoCorrelation_normalize_3chan(ImageGrayd& acorr_1,ImageGrayd& acorr_2,ImageGrayd& acorr_3, double v)
{
	assert (v!=0);
	int T = acorr_1.width();

	for (int y = 0; y < T; ++y)
	{
		for (int x = 0; x < T; ++x)
		{
			acorr_1.pixelAbsolute(x,y) *= v;
			acorr_2.pixelAbsolute(x,y) *= v;
			acorr_3.pixelAbsolute(x,y) *= v;
		}
	}
}


void ASTEX_API shiftDouble (double& d, double s)
{
	// computes -PI < r < PI such that r = d+s mod 2*PI
	// works for -3*PI < d+s < +3*PI
	d += s;
	if (d < -M_PI) d += 2.0*M_PI;
	if (d >  M_PI) d -= 2.0*M_PI;
}


void ASTEX_API phaseShiftFromSpaceShift(ImageSpectrald& phase, double dx, double dy)
{
	// décallage des phases de manière à créer un décallage spatial [dx,dy]
	// hypothèse : la fréquence nulle est en [0,0], i.e. FFTShiftImageFilter n'a pas été appelé
	// hypothèse : phase->GetLargestPossibleRegion().GetIndex() vaut [0,0]

	int NX = phase.width();
	int NY = phase.height();

	for (int j=0; j< NY; ++j)
	{
		for (int i=0; i< NX; ++i)
		{
			double s = - 2.0 * M_PI * fmod ( i*dx/ double(NX) + j*dy/double(NY) , 1.0 );
			assert ( s >= -2*M_PI || s <= 2*M_PI);
			shiftDouble(phase.pixelAbsolute(i,j), s);
		}
	}
}



void ASTEX_API fftForwardModulusAndPhase(const ImageGrayd& input, ImageSpectrald& modulus, ImageSpectrald& phase, bool preserve_energy )
{
	typedef itk::ForwardFFTImageFilter< ImageGrayd::ItkImg, ImageSpectralcd::ItkImg > FFTType;
	FFTType::Pointer fftFilter = FFTType::New();
	fftFilter->SetInput(input.itk());

	typedef itk::FFTShiftImageFilter< ImageSpectralcd::ItkImg, ImageSpectralcd::ItkImg > ShiftFilterType;
	ShiftFilterType::Pointer shiftFilter = ShiftFilterType::New();
	shiftFilter->SetInput(fftFilter->GetOutput());
	shiftFilter->InverseOff();


	typedef itk::ComplexToModulusImageFilter< ImageSpectralcd::ItkImg, ImageSpectrald::ItkImg > ModulusFilterType;
	ModulusFilterType::Pointer modulusFilter = ModulusFilterType::New();
	modulusFilter->SetInput(shiftFilter->GetOutput());

	typedef itk::ShiftScaleImageFilter< ImageSpectrald::ItkImg, ImageSpectrald::ItkImg > RescaleType;
	RescaleType::Pointer rescale = RescaleType::New();
	rescale->SetInput(modulusFilter->GetOutput());
	rescale->SetShift(0.0);

	if(preserve_energy)
	{
		rescale->SetScale(1.0/sqrt(input.width()*input.height()));
	}
	else
	{
		rescale->SetScale(1.0);
	}

	rescale->Update();
	modulus.itk() = rescale->GetOutput();
	modulus.setCenter(modulus.width()/2,modulus.height()/2);

	typedef itk::ComplexToPhaseImageFilter< ImageSpectralcd::ItkImg, ImageSpectrald::ItkImg > PhaseFilterType;
	PhaseFilterType::Pointer phaseFilter = PhaseFilterType::New();
	phaseFilter->SetInput(shiftFilter->GetOutput());

	phaseFilter->Update();

	//	phase.itk() = phaseFilter->GetOutput();
	//	phase.setCenter(phase.width()/2,phase.height()/2);
	phase.update_itk(phaseFilter->GetOutput());
}

void ASTEX_API fftInverseModulusAndPhase(const ImageSpectrald& modulus, const ImageSpectrald& phase, ImageGrayd& output, bool preserve_energy)
{
	typedef itk::MagnitudeAndPhaseToComplexImageFilter< ImageSpectrald::ItkImg, ImageSpectrald::ItkImg, ImageSpectralcd::ItkImg > ToComplexFilterType;
	ToComplexFilterType::Pointer toComplexFilter = ToComplexFilterType::New();
	toComplexFilter->SetInput1(modulus.itk());
	toComplexFilter->SetInput2(phase.itk());

	typedef itk::FFTShiftImageFilter< ImageSpectralcd::ItkImg, ImageSpectralcd::ItkImg > ShiftFilterType;
	ShiftFilterType::Pointer shiftFilter = ShiftFilterType::New();
	shiftFilter->SetInput(toComplexFilter->GetOutput());
	shiftFilter->InverseOn();

	typedef itk::InverseFFTImageFilter< ImageSpectralcd::ItkImg, ImageSpectrald::ItkImg > IFFTType;
	IFFTType::Pointer ifftFilter = IFFTType::New();
	ifftFilter->SetInput(shiftFilter->GetOutput());

	typedef itk::ShiftScaleImageFilter< ImageSpectrald::ItkImg, ImageSpectrald::ItkImg > RescaleType;
	RescaleType::Pointer rescale = RescaleType::New();
	rescale->SetInput(ifftFilter->GetOutput());
	rescale->SetShift(0.0);
	if(preserve_energy)
	{
		rescale->SetScale(sqrt(modulus.width()*modulus.height()));
	}
	else
	{
		rescale->SetScale(1.0);
	}

	rescale->Update();
	//	output.itk() = rescale->GetOutput();
	output.update_itk(rescale->GetOutput());
}

void ASTEX_API welch(const ImageGrayd& input, ImageSpectrald& modulus, unsigned int step)
{
	RegionTraversor w(input.size(),step);
	//	walkInTheImage w (input,step);

	ImageSpectrald modulusAccu(w.getSize()/*,w.getSize()*/,true); // param to true initializes the pixels to 0

	typedef itk::RegionOfInterestImageFilter< ImageGrayd::ItkImg, ImageGrayd::ItkImg > ROIFilterType;
	ROIFilterType::Pointer roiFilter = ROIFilterType::New();

	typedef itk::ForwardFFTImageFilter< ImageGrayd::ItkImg, ImageSpectralcd::ItkImg > FFTType;
	FFTType::Pointer fftFilter = FFTType::New();

	typedef itk::FFTShiftImageFilter< ImageSpectralcd::ItkImg, ImageSpectralcd::ItkImg > ShiftFilterType;
	ShiftFilterType::Pointer shiftFilter = ShiftFilterType::New();

	typedef itk::ComplexToModulusImageFilter< ImageSpectralcd::ItkImg, ImageSpectrald::ItkImg > ModulusFilterType;
	ModulusFilterType::Pointer modulusFilter = ModulusFilterType::New();

	typedef itk::ShiftScaleImageFilter< ImageSpectrald::ItkImg, ImageSpectrald::ItkImg > RescaleType;
	RescaleType::Pointer rescale = RescaleType::New();

	typedef itk::SquareImageFilter<ImageSpectrald::ItkImg, ImageSpectrald::ItkImg> SquareFilterType;
	SquareFilterType::Pointer squareFilter = SquareFilterType::New();

	typedef itk::AddImageFilter<ImageSpectrald::ItkImg> AddFilterType;
	AddFilterType::Pointer addFilter = AddFilterType::New();

	typedef itk::SqrtImageFilter<ImageSpectrald::ItkImg, ImageSpectrald::ItkImg> SqrtFilterType;
	SqrtFilterType::Pointer sqrtFilter = SqrtFilterType::New();


	int count_iterations = 0;
	for (; w.isValid(); ++w)
	{
		count_iterations ++;

		// set the region on [x; x+T-1] * [y; y+T-1]
		roiFilter->SetInput(input.itk());
		roiFilter->SetRegionOfInterest(*w);
		roiFilter->Update();
		roiFilter->GetOutput()->SetOrigin(modulusAccu.itk()->GetOrigin()); // necessary for the addFilter to work

		fftFilter->SetInput(roiFilter->GetOutput());

		shiftFilter->SetInput(fftFilter->GetOutput());
		shiftFilter->InverseOff();

		modulusFilter->SetInput(shiftFilter->GetOutput());

		rescale->SetInput(modulusFilter->GetOutput());
		rescale->SetShift(0.0);
		rescale->SetScale( 1.0 / (double) w.getSize() );

		squareFilter->SetInput(rescale->GetOutput());

		addFilter->SetInput1(modulusAccu.itk());
		addFilter->SetInput2(squareFilter->GetOutput());

		addFilter->Update();
		modulusAccu.itk() = addFilter->GetOutput();
	}

	rescale->SetInput(modulusAccu.itk());
	rescale->SetShift(0.0);
	rescale->SetScale( 1.0 / (double) count_iterations);

	sqrtFilter->SetInput( rescale->GetOutput() );

	sqrtFilter->Update();
	modulus.update_itk(sqrtFilter->GetOutput());
}


void ASTEX_API welch(const ImageGrayd& input, ImageSpectrald& modulus, unsigned int step, int /*sp_size*/)
{
	RegionTraversor w(modulus.size(),step);

	ImageSpectrald modulusAccu(w.getSize()/*,w.getSize()*/,true); // param to true initializes the pixels to 0

	typedef itk::RegionOfInterestImageFilter< ImageGrayd::ItkImg, ImageGrayd::ItkImg > ROIFilterType;
	ROIFilterType::Pointer roiFilter = ROIFilterType::New();

	typedef itk::ForwardFFTImageFilter< ImageGrayd::ItkImg, ImageSpectralcd::ItkImg > FFTType;
	FFTType::Pointer fftFilter = FFTType::New();

	typedef itk::FFTShiftImageFilter< ImageSpectralcd::ItkImg, ImageSpectralcd::ItkImg > ShiftFilterType;
	ShiftFilterType::Pointer shiftFilter = ShiftFilterType::New();

	typedef itk::ComplexToModulusImageFilter< ImageSpectralcd::ItkImg, ImageSpectrald::ItkImg > ModulusFilterType;
	ModulusFilterType::Pointer modulusFilter = ModulusFilterType::New();

	typedef itk::ShiftScaleImageFilter< ImageSpectrald::ItkImg, ImageSpectrald::ItkImg > RescaleType;
	RescaleType::Pointer rescale = RescaleType::New();

	typedef itk::SquareImageFilter<ImageSpectrald::ItkImg, ImageSpectrald::ItkImg> SquareFilterType;
	SquareFilterType::Pointer squareFilter = SquareFilterType::New();

	typedef itk::AddImageFilter<ImageSpectrald::ItkImg> AddFilterType;
	AddFilterType::Pointer addFilter = AddFilterType::New();

	typedef itk::SqrtImageFilter<ImageSpectrald::ItkImg, ImageSpectrald::ItkImg> SqrtFilterType;
	SqrtFilterType::Pointer sqrtFilter = SqrtFilterType::New();


	int count_iterations = 0;
	for (; w.isValid(); ++w)
	{
		count_iterations ++;

		// set the region on [x; x+T-1] * [y; y+T-1]
		roiFilter->SetInput(input.itk());
		roiFilter->SetRegionOfInterest(*w);
		roiFilter->Update();
		roiFilter->GetOutput()->SetOrigin(modulusAccu.itk()->GetOrigin()); // necessary for the addFilter to work

		fftFilter->SetInput(roiFilter->GetOutput());

		shiftFilter->SetInput(fftFilter->GetOutput());
		shiftFilter->InverseOff();

		modulusFilter->SetInput(shiftFilter->GetOutput());

		rescale->SetInput(modulusFilter->GetOutput());
		rescale->SetShift(0.0);
		rescale->SetScale( 1.0 / (double) w.getSize() );

		squareFilter->SetInput(rescale->GetOutput());

		addFilter->SetInput1(modulusAccu.itk());
		addFilter->SetInput2(squareFilter->GetOutput());

		addFilter->Update();
		modulusAccu.itk() = addFilter->GetOutput();
	}

	rescale->SetInput(modulusAccu.itk());
	rescale->SetShift(0.0);
	rescale->SetScale( 1.0 / (double) count_iterations);

	sqrtFilter->SetInput( rescale->GetOutput() );

	sqrtFilter->Update();
	modulus.update_itk(sqrtFilter->GetOutput());
}







void ASTEX_API setPower(ImageGrayd& image, double power)
{

	double s = std::sqrt(power/getPower(image));

	image.for_all_pixels([&s](double& p)
	{
		p *= s;
	});
}


void ASTEX_API down_sampling(const ImageSpectrald& modulus, ImageSpectrald& output)
{
	assert(modulus.width() == modulus.height());
	int T = modulus.width();
	assert(output.width() == T/2 && output.height() == T/2);

	ImageSpectrald psd(T);

	for (int y=0; y<T; y++)
	{
		for (int x = 0; x< T; x++)
		{
			const double v = modulus.pixelAbsolute(x,y);
			psd.pixelAbsolute(x,y) = v*v;
		}
	}

	for (int y=1; y<T; y+=2)
	{
		for (int x = 1; x< T; x+=2)
		{
			psd.pixelAbsolute(x,y) += psd.pixelAbsolute((x-1)% T,(y)%T) / 2.0;
			psd.pixelAbsolute(x,y) += psd.pixelAbsolute((x+1)% T,(y)%T) / 2.0;
			psd.pixelAbsolute(x,y) += psd.pixelAbsolute((x)% T,(y-1)%T) / 2.0;
			psd.pixelAbsolute(x,y) += psd.pixelAbsolute((x)% T,(y+1)%T) / 2.0;

			psd.pixelAbsolute(x,y) += psd.pixelAbsolute((x-1)% T,(y-1)%T) / 4.0;
			psd.pixelAbsolute(x,y) += psd.pixelAbsolute((x+1)% T,(y-1)%T) / 4.0;
			psd.pixelAbsolute(x,y) += psd.pixelAbsolute((x-1)% T,(y+1)%T) / 4.0;
			psd.pixelAbsolute(x,y) += psd.pixelAbsolute((x+1)% T,(y+1)%T) / 4.0;
		}
	}

	for (int y=0; y<T/2; y++)
	{
		for (int x=0; x< T/2; x++)
		{
			output.pixelAbsolute(x,y) = std::sqrt(psd.pixelAbsolute(2*x+1,2*y+1)) / 2.0;
		}
	}
}


void ASTEX_API phaseOfRegion(const ImageGrayd& input, const Region& region, ImageSpectrald& phase )
{
	typedef itk::RegionOfInterestImageFilter< ImageGrayd::ItkImg, ImageGrayd::ItkImg > ROIFilterType;
	ROIFilterType::Pointer roiFilter = ROIFilterType::New();
	roiFilter->SetInput(input.itk());
	roiFilter->SetRegionOfInterest(region);
	roiFilter->Update();

	typedef itk::ForwardFFTImageFilter< ImageGrayd::ItkImg, ImageSpectralcd::ItkImg > FFTType;
	FFTType::Pointer fftFilter = FFTType::New();
	fftFilter->SetInput(roiFilter->GetOutput());

	typedef itk::FFTShiftImageFilter< ImageSpectralcd::ItkImg, ImageSpectralcd::ItkImg > ShiftFilterType;
	ShiftFilterType::Pointer shiftFilter = ShiftFilterType::New();
	shiftFilter->SetInput(fftFilter->GetOutput());
	shiftFilter->InverseOff();

	typedef itk::ComplexToPhaseImageFilter< ImageSpectralcd::ItkImg, ImageSpectrald::ItkImg > PhaseFilterType;
	PhaseFilterType::Pointer phaseFilter = PhaseFilterType::New();
	phaseFilter->SetInput(shiftFilter->GetOutput());
	phaseFilter->Update();
	phase.itk() = phaseFilter->GetOutput();
}

void ASTEX_API phaseAverage_UsingNormalizedComplexAverage(const ImageGrayd& input, ImageSpectrald& phase, bool activate, unsigned int step)
{
	RegionTraversor w (input.size(),step);

	ImageGraycd accu(w.getSize(),w.getSize(),true); // param to true initializes the pixels to 0

	ImageSpectrald constantModulus(w.getSize()/*,w.getSize()*/,true); // param to true initializes the pixels to 0
	constantModulus.itk()->FillBuffer(1);

	typedef itk::RegionOfInterestImageFilter< ImageGrayd::ItkImg, ImageGrayd::ItkImg > ROIFilterType;
	ROIFilterType::Pointer roiFilter = ROIFilterType::New();

	typedef itk::CyclicShiftImageFilter< ImageGrayd::ItkImg, ImageGrayd::ItkImg > CyclicShiftFilterType;
	CyclicShiftFilterType::Pointer cyclicShiftFilter = CyclicShiftFilterType::New();

	typedef itk::ForwardFFTImageFilter< ImageGrayd::ItkImg, ImageSpectralcd::ItkImg > FFTType;
	FFTType::Pointer fftFilter = FFTType::New();

	typedef itk::FFTShiftImageFilter< ImageSpectralcd::ItkImg, ImageSpectralcd::ItkImg > ShiftFilterType;
	ShiftFilterType::Pointer fftShiftFilter = ShiftFilterType::New();

	typedef itk::ComplexToPhaseImageFilter< ImageSpectralcd::ItkImg, ImageSpectrald::ItkImg > PhaseFilterType;
	PhaseFilterType::Pointer phaseFilter = PhaseFilterType::New();

	typedef itk::MagnitudeAndPhaseToComplexImageFilter< ImageSpectrald::ItkImg, ImageSpectrald::ItkImg, ImageSpectralcd::ItkImg > ToComplexFilterType;
	ToComplexFilterType::Pointer toComplexFilter = ToComplexFilterType::New();

	typedef itk::AddImageFilter<ImageSpectralcd::ItkImg> AddFilterType;
	AddFilterType::Pointer addFilter = AddFilterType::New();


	int count_iterations = 0;
	for (; w.isValid(); ++w)
	{
		count_iterations ++;

		roiFilter->SetInput(input.itk());
		roiFilter->SetRegionOfInterest(*w);
		roiFilter->Update();
		roiFilter->GetOutput()->SetOrigin(accu.itk()->GetOrigin()); // necessary for the addFilter to work

		cyclicShiftFilter->SetInput(roiFilter->GetOutput());
		if (activate)
		{
			Offset offset;
			offset[0]=-w.getX();
			offset[1]=-w.getY();
			cyclicShiftFilter->SetShift(offset);
		}
		else
		{
			/*Offset offset({{0,0}});*/
			Offset offset;
			offset[0]=0;
			offset[1]=0;
			cyclicShiftFilter->SetShift(offset);
		}

		fftFilter->SetInput(cyclicShiftFilter->GetOutput());

		fftShiftFilter->SetInput(fftFilter->GetOutput());
		fftShiftFilter->InverseOff();

		phaseFilter->SetInput(fftShiftFilter->GetOutput());

		toComplexFilter->SetInput1(constantModulus.itk());
		toComplexFilter->SetInput2(phaseFilter->GetOutput());

		addFilter->SetInput1(accu.itk());
		addFilter->SetInput2(toComplexFilter->GetOutput());

		addFilter->Update();
		accu.itk() = addFilter->GetOutput();
	}

	phaseFilter->SetInput(accu.itk());


	phaseFilter->Update();
	phase.update_itk(phaseFilter->GetOutput());
}



void ASTEX_API phaseAverage_Naive(const ImageGrayd& input, ImageSpectrald& phase, bool activate, unsigned int step)
{
	RegionTraversor w (input.size(),step);

	ImageSpectrald phaseAccu(w.getSize()/*,w.getSize()*/,true); // param to true initializes the pixels to 0

	typedef itk::RegionOfInterestImageFilter< ImageGrayd::ItkImg, ImageGrayd::ItkImg > ROIFilterType;
	ROIFilterType::Pointer roiFilter = ROIFilterType::New();

	typedef itk::CyclicShiftImageFilter< ImageGrayd::ItkImg, ImageGrayd::ItkImg > CyclicShiftFilterType;
	CyclicShiftFilterType::Pointer cyclicShiftFilter = CyclicShiftFilterType::New();

	typedef itk::ForwardFFTImageFilter< ImageGrayd::ItkImg, ImageSpectralcd::ItkImg > FFTType;
	FFTType::Pointer fftFilter = FFTType::New();

	typedef itk::FFTShiftImageFilter< ImageSpectrald::ItkImg, ImageSpectrald::ItkImg > ShiftFilterType;
	ShiftFilterType::Pointer shiftFilter = ShiftFilterType::New();

	typedef itk::ComplexToPhaseImageFilter< ImageSpectralcd::ItkImg, ImageSpectrald::ItkImg > PhaseFilterType;
	PhaseFilterType::Pointer phaseFilter = PhaseFilterType::New();

	typedef itk::AddImageFilter<ImageSpectrald::ItkImg> AddFilterType;
	AddFilterType::Pointer addFilter = AddFilterType::New();

	int count_iterations = 0;
	for (; w.isValid(); ++w)
	{
		count_iterations ++;

		roiFilter->SetInput(input.itk());
		roiFilter->SetRegionOfInterest(*w);
		roiFilter->Update();
		roiFilter->GetOutput()->SetOrigin(phaseAccu.itk()->GetOrigin()); // necessary for the addFilter to work

		cyclicShiftFilter->SetInput(roiFilter->GetOutput());
		if (activate)
		{
			Offset offset/*({{-w.getX(),-w.getY()}})*/;
			offset[0]=-w.getX();
			offset[1]=-w.getY();
			cyclicShiftFilter->SetShift(offset);
		}
		else
		{
			Offset offset/*({{0,0}})*/;
			offset[0]=0;
			offset[1]=0;
			cyclicShiftFilter->SetShift(offset);
		}

		fftFilter->SetInput(cyclicShiftFilter->GetOutput());

		phaseFilter->SetInput(fftFilter->GetOutput());

		shiftFilter->SetInput(phaseFilter->GetOutput());
		shiftFilter->InverseOff();

		addFilter->SetInput1(phaseAccu.itk());
		addFilter->SetInput2(shiftFilter->GetOutput());

		addFilter->Update();
		phaseAccu.itk() = addFilter->GetOutput();
	}

	typedef itk::ShiftScaleImageFilter< ImageSpectrald::ItkImg, ImageSpectrald::ItkImg > RescaleType;
	RescaleType::Pointer rescale = RescaleType::New();
	rescale->SetInput(phaseAccu.itk());
	rescale->SetShift(0.0);
	rescale->SetScale( 1.0 / (double) count_iterations);

	rescale->Update();
	phase.update_itk(rescale->GetOutput());
}


void ASTEX_API phaseAverageUsingPhaseShift(const ImageGrayd& input, ImageSpectrald& phase, bool activate, unsigned int step)
{
	RegionTraversor w (input.size(),step);

	ImageSpectrald phaseAccu(w.getSize()/*,w.getSize()*/,true); // param to true initializes the pixels to 0

	typedef itk::RegionOfInterestImageFilter< ImageGrayd::ItkImg, ImageGrayd::ItkImg > ROIFilterType;
	ROIFilterType::Pointer roiFilter = ROIFilterType::New();

	typedef itk::ForwardFFTImageFilter< ImageGrayd::ItkImg, ImageSpectralcd::ItkImg > FFTType;
	FFTType::Pointer fftFilter = FFTType::New();

	typedef itk::FFTShiftImageFilter< ImageSpectrald::ItkImg, ImageSpectrald::ItkImg > ShiftFilterType;
	ShiftFilterType::Pointer shiftFilter = ShiftFilterType::New();

	typedef itk::ComplexToPhaseImageFilter< ImageSpectralcd::ItkImg, ImageSpectrald::ItkImg > PhaseFilterType;
	PhaseFilterType::Pointer phaseFilter = PhaseFilterType::New();

	typedef itk::AddImageFilter<ImageSpectrald::ItkImg> AddFilterType;
	AddFilterType::Pointer addFilter = AddFilterType::New();

	int count_iterations = 0;
	for (; w.isValid(); ++w)
	{
		count_iterations ++;

		roiFilter->SetInput(input.itk());
		roiFilter->SetRegionOfInterest(*w);
		roiFilter->Update();
		roiFilter->GetOutput()->SetOrigin(phaseAccu.itk()->GetOrigin()); // necessary for the addFilter to work

		fftFilter->SetInput(roiFilter->GetOutput());

		phaseFilter->SetInput(fftFilter->GetOutput());

		if (activate)
		{
			phaseFilter->Update();
			ImageSpectrald im (phaseFilter->GetOutput());
			phaseShiftFromSpaceShift(im,-w.getX(),-w.getY()); // must be applied before shiftFilter

		}

		shiftFilter->SetInput(phaseFilter->GetOutput());
		shiftFilter->InverseOff();

		addFilter->SetInput1(phaseAccu.itk());
		addFilter->SetInput2(shiftFilter->GetOutput());

		addFilter->Update();
		phaseAccu.itk() = addFilter->GetOutput();
	}

	typedef itk::ShiftScaleImageFilter< ImageSpectrald::ItkImg, ImageSpectrald::ItkImg > RescaleType;
	RescaleType::Pointer rescale = RescaleType::New();
	rescale->SetInput(phaseAccu.itk());
	rescale->SetShift(0.0);
	rescale->SetScale( 1.0 / (double) count_iterations);

	rescale->Update();
	phase.update_itk(rescale->GetOutput());
}


void ASTEX_API RPnoise_mosaic(const ImageSpectrald& modulus, ImageGrayd& result, int blending_size,bool call_srand)
{
	const int t_w = modulus.width();
	const int t_h = modulus.height();
	assert(t_w == t_h);
	assert(blending_size < t_w / 2 && blending_size < t_h / 2 && blending_size >= 0);

	const int im_w = result.width();
	const int im_h = result.height();
	for (int x=0; x< im_w; ++x) for (int y=0; y<im_h; ++y) result.pixelAbsolute(x,y) = 0.0;

	if (call_srand) srand (time(NULL));

	ImageGrayd blending(t_w,t_h,true);
	for (int x=0; x< t_w; ++x)
	{
		double vx = 1.0;
		if (x < blending_size)
			vx = double(x+1) / double(blending_size+1);
		if (x >= t_w - blending_size)
			vx = double(t_w - x) / double(blending_size+1);

		for (int y=0; y<t_h; ++y)
		{
			double vy = 1.0;
			if (y < blending_size)
				vy = double(y+1) / double(blending_size+1);
			if (y >= t_h - blending_size)
				vy = double(t_h - y) / double(blending_size+1);
			blending.pixelAbsolute(x,y) = vx * vy;
		}
	}



	for (int y=-blending_size; y - blending_size < im_h; y+=t_h-blending_size)
	{
		for (int x=-blending_size; x - blending_size < im_w; x+=t_w-blending_size)
		{ // compute + copy + blend the t_w*t_h block starting at (x,y)

			ImageSpectrald phase(t_w,/*t_h,*/true);
			randomPhase(phase, [](int,int){return true;}, false);

			ImageGrayd local_noise;
			fftInverseModulusAndPhase(modulus, phase, local_noise);

			for (int j=0; j<t_h; j++)
			{
				for (int i = 0; i<t_w; ++i)
				{
					const int pi = x+i, pj = y+j;
					if (0 <= pi && pi < im_w && 0 <= pj && pj < im_h)
						result.pixelAbsolute(pi,pj) += blending.pixelAbsolute(i,j) * local_noise.pixelAbsolute(i,j);
				}
			}
		}
	}

}



void ASTEX_API RPnoise_mosaic_periodique(const ImageSpectrald& modulus, ImageGrayd& result)
{
	const int t_w = modulus.width();
	const int t_h = modulus.height();
	assert(t_w == t_h);

	const int im_w = result.width();
	const int im_h = result.height();
	for (int x=0; x< im_w; ++x) for (int y=0; y<im_h; ++y) result.pixelAbsolute(x,y) = 0.0;

	ImageSpectrald phase(t_w,/*t_h,*/true);
	randomPhase(phase, [](int,int){return true;}, false);

	ImageGrayd local_noise;
	fftInverseModulusAndPhase(modulus, phase, local_noise);

	for (int y=0; y < im_h; y+=t_h )
	{
		for (int x=0 ; x < im_w; x+=t_w )
		{
			for (int j=0; j<t_h; j++)
			{
				for (int i = 0; i<t_w; ++i)
				{
					const int pi = x+i, pj = y+j;
					if (0 <= pi && pi < im_w && 0 <= pj && pj < im_h)
						result.pixelAbsolute(pi,pj) = local_noise.pixelAbsolute(i,j);
				}
			}
		}
	}

}

void ASTEX_API RPnoise_mosaic_bandes(const ImageSpectrald& modulus, ImageGrayd& result)
{

	const int t_w = modulus.width();
	const int t_h = modulus.height();
	assert(t_w == t_h);

	const int im_w = result.width();
	const int im_h = result.height();
	for (int x=0; x< im_w; ++x) for (int y=0; y<im_h; ++y) result.pixelAbsolute(x,y) = 0.0;

	ImageGrayd result_0(im_w,im_h,true);
	ImageGrayd result_1(im_w,im_h,true);
	ImageGrayd result_2(im_w,im_h,true);


	ImageSpectrald modulus_0(t_w/*,t_h*/,true);
	ImageSpectrald modulus_1(t_w/*,t_h*/,true);
	ImageSpectrald modulus_2(t_w/*,t_h*/,true);

	double min,max;
	min = modulus.pixelRelative(0,0);
	max = modulus.pixelRelative(0,0);

	int seuil_01 = 0;
	int seuil_12 = 0;

	std::vector<double> e(t_h/2+1,0.0);
	//   e.resize(t_h/2+1,0);
	double esum  = 0;

	for (int j=-t_h/2; j<t_h/2; j++)
	{
		for (int i = -t_w/2; i<t_w/2; ++i)
		{
			const double v = modulus.pixelRelative(i,j);
			e[std::max(abs(i) , abs(j))] += v*v;
			esum += v*v;
		}
	}
	//   std::cout << "SUM" << esum << std::endl;

	esum -= e[0]; // do not take into account the null frequency
	double ecum  = 0; // e[0];
	int j=0;
	while (ecum < esum/3)
	{
		++j;
		ecum += e[j];
	}
	seuil_12 = j;
	while (ecum < esum*2/3)
	{
		++j;
		ecum += e[j];
	}
	seuil_01 = j;


	for (int j=-t_h/2; j<t_h/2; j++)
	{
		for (int i = -t_w/2; i<t_w/2; ++i)
		{
			const double v = modulus.pixelRelative(i,j);
			max = std::max(max,v);
			min = std::min(min,v);

			if( std::max(abs(i) , abs(j)) > seuil_01)
			{
				modulus_0.pixelRelative(i,j) = v;
			}
			else if(std::max(abs(i) , abs(j)) > seuil_12)
			{
				modulus_1.pixelRelative(i,j) = v;
			}
			else
			{
				modulus_2.pixelRelative(i,j) = v;
			}
		}
	}
	//    std::cout << "ok" << (seuil_01 + t_h / 2) << " " << seuil_01 << " " <<  seuil_12 << std::endl;

	int blend_0 = t_h*2 / (seuil_01 + t_h / 2); // mid stratum
	int blend_1 = t_h / seuil_01; // max frequency = min period
	int blend_2 = t_h / seuil_12; // max frequency = min period

	blend_0 = std::min(blend_0, t_h/4); // for safety reasons
	blend_1 = std::min(blend_1, t_h/4); // for safety reasons
	blend_2 = std::min(blend_2, t_h/4); // for safety reasons

	RPnoise_mosaic(modulus_0,result_0,blend_0);
	RPnoise_mosaic(modulus_1,result_1,blend_1);
	RPnoise_mosaic(modulus_2,result_2,blend_2);

	for (int y=0; y < im_h; y++ )
	{
		for (int x=0 ; x < im_w; x++ )
		{
			result.pixelAbsolute(x,y) = result_0.pixelAbsolute(x,y) + result_1.pixelAbsolute(x,y) + result_2.pixelAbsolute(x,y) ;
		}
	}
}




double ASTEX_API distance_spectrum_to_spectrum_linear_weights(const ImageSpectrald& sp1, const ImageSpectrald& sp2)
{
	return std::sqrt(distance_squared_spectrum_to_spectrum_linear_weights(sp1,sp2) );
}

double ASTEX_API distance_spectrum_to_spectrum_uniform_weights(const ImageSpectrald& sp1, const ImageSpectrald& sp2)
{
	return std::sqrt(distance_squared_spectrum_to_spectrum_uniform_weights(sp1,sp2) );
}

double ASTEX_API distance_squared_spectrum_to_spectrum_linear_weights(const ImageSpectrald& sp1, const ImageSpectrald& sp2)
{
	const int T = sp1.width();
	assert ( sp1.height()==T && sp2.width()==T && sp2.height()==T );
	assert ( sp1.getCenter()[0]==T/2 && sp1.getCenter()[1]==T/2 && sp2.getCenter()[0]==T/2 && sp2.getCenter()[1]==T/2);

	double r = 0.0;
	for_indices(-T/2, T/2, -T/2, T/2, [&] (int x, int y)
	{
		const double weight = std::max( 1.0, std::sqrt(x*x+y*y) ) / double(T);
		const double ds = sp1.pixelRelative(x,y) - sp2.pixelRelative(x,y);
		r += ds * ds * weight;
	});
	return r;
}

double ASTEX_API distance_squared_spectrum_to_spectrum_uniform_weights(const ImageSpectrald& sp1, const ImageSpectrald& sp2)
{
	const int T = sp1.width();
	assert ( sp1.height()==T && sp2.width()==T && sp2.height()==T );
	assert ( sp1.getCenter()[0]==T/2 && sp1.getCenter()[1]==T/2 && sp2.getCenter()[0]==T/2 && sp2.getCenter()[1]==T/2);

	double r = 0.0;
	for_indices(-T/2, T/2, -T/2, T/2, [&] (int x, int y)
	{
		const double weight = 1.0 / double(T);
		const double ds = sp1.pixelRelative(x,y) - sp2.pixelRelative(x,y);
		r += ds * ds * weight;
	});

	return r;
}


double ASTEX_API dot_product_spectrum_linear_weights(const ImageSpectrald& sp1, const ImageSpectrald& sp2)
{
	const int T = sp1.width();
	assert ( sp1.height()==T && sp2.width()==T && sp2.height()==T );
	assert ( sp1.getCenter()[0]==T/2 && sp1.getCenter()[1]==T/2 && sp2.getCenter()[0]==T/2 && sp2.getCenter()[1]==T/2);

	double r = 0.0;
	for_indices(-T/2, T/2, -T/2, T/2, [&] (int x, int y)
	{
		const double weight = std::max( 1.0, std::sqrt(x*x+y*y) ) / double(T);
		r += sp1.pixelRelative(x,y) * sp2.pixelRelative(x,y) * weight;
	});

	return r;
}

double ASTEX_API dot_product_spectrum_uniform_weights(const ImageSpectrald& sp1, const ImageSpectrald& sp2)
{
	const int T = sp1.width();
	assert ( sp1.height()==T && sp2.width()==T && sp2.height()==T );
	assert ( sp1.getCenter()[0]==T/2 && sp1.getCenter()[1]==T/2 && sp2.getCenter()[0]==T/2 && sp2.getCenter()[1]==T/2);

	double r = 0.0;
	for_indices(-T/2, T/2, -T/2, T/2, [&] (int x, int y)
	{
		const double weight = 1.0 / double(T);
		r += sp1.pixelRelative(x,y) * sp2.pixelRelative(x,y) * weight;
	});

	return r;
}

void spectrum_projector_2D::add (const ImageSpectrald& sp)
{
	m_S.push_back(sp);
}

void spectrum_projector_2D::precompute_system()
{
	assert(m_S.size() == 2);

	m_A[0][0] = dot_product_spectrum_linear_weights(m_S[0],m_S[0]);
	m_A[0][1] = dot_product_spectrum_linear_weights(m_S[0],m_S[1]);
	m_A[1][0] = m_A[0][1]; // symmetry of dot product
	m_A[1][1] = dot_product_spectrum_linear_weights(m_S[1],m_S[1]);

	double det = m_A[0][0]*m_A[1][1] - m_A[1][0]*m_A[0][1];
	m_A_inverse[0][0] = m_A[1][1] / det;
	m_A_inverse[1][1] = m_A[0][0] / det;
	m_A_inverse[0][1] = - m_A[0][1] / det;
	m_A_inverse[1][0] = - m_A[1][0] / det;
}

void spectrum_projector_2D::project (const ImageSpectrald& sp, std::vector<double>& coord)
{
	coord.resize(2);
	double b [2];
	b[0] = dot_product_spectrum_linear_weights(sp,m_S[0]);
	b[1] = dot_product_spectrum_linear_weights(sp,m_S[1]);
	coord[0] = m_A_inverse[0][0]*b[0] + m_A_inverse[0][1]*b[1];
	coord[1] = m_A_inverse[1][0]*b[0] + m_A_inverse[1][1]*b[1];
}

void spectrum_projector_2D::precompute_system_uniform_weights()
{
	assert(m_S.size() == 2);

	m_A[0][0] = dot_product_spectrum_uniform_weights(m_S[0],m_S[0]);
	m_A[0][1] = dot_product_spectrum_uniform_weights(m_S[0],m_S[1]);
	m_A[1][0] = m_A[0][1]; // symmetry of dot product
	m_A[1][1] = dot_product_spectrum_uniform_weights(m_S[1],m_S[1]);

	double det = m_A[0][0]*m_A[1][1] - m_A[1][0]*m_A[0][1];
	m_A_inverse[0][0] = m_A[1][1] / det;
	m_A_inverse[1][1] = m_A[0][0] / det;
	m_A_inverse[0][1] = - m_A[0][1] / det;
	m_A_inverse[1][0] = - m_A[1][0] / det;
}

void spectrum_projector_2D::project_uniform_weights (const ImageSpectrald& sp, std::vector<double>& coord)
{
	coord.resize(2);
	double b [2];
	b[0] = dot_product_spectrum_uniform_weights(sp,m_S[0]);
	b[1] = dot_product_spectrum_uniform_weights(sp,m_S[1]);
	coord[0] = m_A_inverse[0][0]*b[0] + m_A_inverse[0][1]*b[1];
	coord[1] = m_A_inverse[1][0]*b[0] + m_A_inverse[1][1]*b[1];
}


void spectrum_projector_2D::reproject_coordinates_to_positive (std::vector<double>& coord)
{
	if (coord[0]<0)
	{
		reproject(coord,1);
		if (coord[1]<0) coord[1]=0;
	}
	else if (coord[1]<0)
	{
		reproject(coord,0);
		if (coord[0]<0) coord[0]=0;
	}
}

void spectrum_projector_2D::reproject_coordinates_to_0_1 (std::vector<double>& coord)
{
	assert(coord[0]>=0 && coord[1]>=0);
	if (coord[0]>1 || coord[1]>1)
	{
		const double sum = coord[0] + coord[1];
		coord[0] /= sum;
		coord[1] /= sum;
	}
}

void spectrum_projector_2D::project_with_positive_coordinates(const ImageSpectrald& sp, std::vector<double>& coord)
{
	project(sp,coord);
	reproject_coordinates_to_positive(coord);
}

double spectrum_projector_2D::dot_product (std::vector<double>& c1 , std::vector<double>& c2)
{
	assert(c1.size()==2 && c2.size()==2);
	return m_A[0][0]*c1[0]*c2[0] + m_A[0][1]*c1[0]*c2[1] + m_A[1][0]*c1[1]*c2[0]  + m_A[1][1]*c1[1]*c2[1];
}

double spectrum_projector_2D::distance_spectrum_to_spectrum (std::vector<double>& c1 , std::vector<double>& c2)
{
	std::vector<double> d;

	d.push_back(c1[0]-c2[0]);
	d.push_back(c1[1]-c2[1]);

	return std::sqrt(this->dot_product(d,d));
}

void spectrum_projector_2D::reproject (std::vector<double>& coord, const unsigned int x)
{
	const unsigned int y = 1-x;
	coord[x] += coord[y] * m_A[x][y] / m_A[x][x];
	coord[y] = 0;
}

void spectrum_projector_3D::add (const ImageSpectrald& sp)
{
	m_S.push_back(sp);
}

void spectrum_projector_3D::precompute_system()
{
	assert(m_S.size() == 3);

	m_A[0][0] = dot_product_spectrum_linear_weights(m_S[0],m_S[0]);
	m_A[0][1] = dot_product_spectrum_linear_weights(m_S[0],m_S[1]);
	m_A[1][0] = m_A[0][1]; // symmetry of dot product
	m_A[0][2] = dot_product_spectrum_linear_weights(m_S[0],m_S[2]);
	m_A[2][0] = m_A[0][2]; // symmetry of dot product
	m_A[1][1] = dot_product_spectrum_linear_weights(m_S[1],m_S[1]);
	m_A[1][2] = dot_product_spectrum_linear_weights(m_S[1],m_S[2]);
	m_A[2][1] = m_A[1][2]; // symmetry of dot product
	m_A[2][2] = dot_product_spectrum_linear_weights(m_S[2],m_S[2]);

	double det = m_A[0][0]*m_A[1][1]*m_A[2][2] + m_A[0][1]*m_A[1][2]*m_A[2][0] + m_A[0][2]*m_A[1][0]*m_A[2][1] - m_A[0][0]*m_A[1][2]*m_A[2][1] - m_A[0][1]*m_A[1][0]*m_A[2][2] - m_A[0][2]*m_A[1][1]*m_A[2][0];
	m_A_inverse[0][0] = (m_A[1][1] * m_A[2][2] - m_A[2][1] * m_A[1][2]) / det;
	m_A_inverse[0][1] = - (m_A[1][0] * m_A[2][2] - m_A[2][0] * m_A[1][2]) / det;
	m_A_inverse[1][0] = m_A_inverse[0][1]; // symmetry
	m_A_inverse[0][2] = (m_A[1][0] * m_A[2][1] - m_A[2][0] * m_A[1][1]) / det;
	m_A_inverse[2][0] = m_A_inverse[0][2]; // symmetry
	m_A_inverse[1][1] = (m_A[0][0] * m_A[2][2] - m_A[2][0] * m_A[0][2]) / det;
	m_A_inverse[1][2] = - (m_A[0][0] * m_A[2][1] - m_A[2][0] * m_A[0][1]) / det;
	m_A_inverse[2][1] = m_A_inverse[1][2]; // symmetry
	m_A_inverse[2][2] = (m_A[0][0] * m_A[1][1] - m_A[1][0] * m_A[0][1]) / det;
}

void spectrum_projector_3D::project (const ImageSpectrald& sp, std::vector<double>& coord)
{
	coord.resize(3);
	double b [3];
	b[0] = dot_product_spectrum_linear_weights(sp,m_S[0]);
	b[1] = dot_product_spectrum_linear_weights(sp,m_S[1]);
	b[2] = dot_product_spectrum_linear_weights(sp,m_S[2]);
	coord[0] = m_A_inverse[0][0]*b[0] + m_A_inverse[0][1]*b[1] + m_A_inverse[0][2]*b[2];
	coord[1] = m_A_inverse[1][0]*b[0] + m_A_inverse[1][1]*b[1] + m_A_inverse[1][2]*b[2];
	coord[2] = m_A_inverse[2][0]*b[0] + m_A_inverse[2][1]*b[1] + m_A_inverse[2][2]*b[2];
}

void spectrum_projector_3D::precompute_system_uniform_weights()
{
	assert(m_S.size() == 3);

	m_A[0][0] = dot_product_spectrum_uniform_weights(m_S[0],m_S[0]);
	m_A[0][1] = dot_product_spectrum_uniform_weights(m_S[0],m_S[1]);
	m_A[1][0] = m_A[0][1]; // symmetry of dot product
	m_A[0][2] = dot_product_spectrum_uniform_weights(m_S[0],m_S[2]);
	m_A[2][0] = m_A[0][2]; // symmetry of dot product
	m_A[1][1] = dot_product_spectrum_uniform_weights(m_S[1],m_S[1]);
	m_A[1][2] = dot_product_spectrum_uniform_weights(m_S[1],m_S[2]);
	m_A[2][1] = m_A[1][2]; // symmetry of dot product
	m_A[2][2] = dot_product_spectrum_uniform_weights(m_S[2],m_S[2]);

	double det = m_A[0][0]*m_A[1][1]*m_A[2][2] + m_A[0][1]*m_A[1][2]*m_A[2][0] + m_A[0][2]*m_A[1][0]*m_A[2][1] - m_A[0][0]*m_A[1][2]*m_A[2][1] - m_A[0][1]*m_A[1][0]*m_A[2][2] - m_A[0][2]*m_A[1][1]*m_A[2][0];
	m_A_inverse[0][0] = (m_A[1][1] * m_A[2][2] - m_A[2][1] * m_A[1][2]) / det;
	m_A_inverse[0][1] = - (m_A[1][0] * m_A[2][2] - m_A[2][0] * m_A[1][2]) / det;
	m_A_inverse[1][0] = m_A_inverse[0][1]; // symmetry
	m_A_inverse[0][2] = (m_A[1][0] * m_A[2][1] - m_A[2][0] * m_A[1][1]) / det;
	m_A_inverse[2][0] = m_A_inverse[0][2]; // symmetry
	m_A_inverse[1][1] = (m_A[0][0] * m_A[2][2] - m_A[2][0] * m_A[0][2]) / det;
	m_A_inverse[1][2] = - (m_A[0][0] * m_A[2][1] - m_A[2][0] * m_A[0][1]) / det;
	m_A_inverse[2][1] = m_A_inverse[1][2]; // symmetry
	m_A_inverse[2][2] = (m_A[0][0] * m_A[1][1] - m_A[1][0] * m_A[0][1]) / det;
}

void spectrum_projector_3D::project_uniform_weights (const ImageSpectrald& sp, std::vector<double>& coord)
{
	coord.resize(3);
	double b [3];
	b[0] = dot_product_spectrum_uniform_weights(sp,m_S[0]);
	b[1] = dot_product_spectrum_uniform_weights(sp,m_S[1]);
	b[2] = dot_product_spectrum_uniform_weights(sp,m_S[2]);
	coord[0] = m_A_inverse[0][0]*b[0] + m_A_inverse[0][1]*b[1] + m_A_inverse[0][2]*b[2];
	coord[1] = m_A_inverse[1][0]*b[0] + m_A_inverse[1][1]*b[1] + m_A_inverse[1][2]*b[2];
	coord[2] = m_A_inverse[2][0]*b[0] + m_A_inverse[2][1]*b[1] + m_A_inverse[2][2]*b[2];
}

void spectrum_projector_3D::reproject_coordinates_to_positive (std::vector<double>& coord)
{
	if (coord[0]<0)
	{
		if (coord[1]<0)
		{
			if (coord[2]<0)
				; // do nothing, all coordinates will be set to zero
			else
				reproject(coord,2);
		}
		else
		{
			if (coord[2]<0)
				reproject(coord,1);
			else
				reproject(coord,1,2);
		}

	}
	else if (coord[1]<0)
	{
		if (coord[2]<0)
			reproject(coord,0);
		else
			reproject(coord,0,2);
	}
	else if (coord[2]<0)
	{
		reproject(coord,0,1);
	}

	if (coord[0]<0) coord[0]=0;
	if (coord[1]<0) coord[1]=0;
	if (coord[2]<0) coord[2]=0;

}

void spectrum_projector_3D::reproject_coordinates_to_0_1 (std::vector<double>& coord)
{
	assert(coord[0]>=0 && coord[1]>=0 && coord[2]>=0);
	if (coord[0]>1 || coord[1]>1 || coord[2]>1)
	{
		const double sum = coord[0] + coord[1] + coord[2];
		coord[0] /= sum;
		coord[1] /= sum;
		coord[2] /= sum;
	}
}

void spectrum_projector_3D::project_with_positive_coordinates(const ImageSpectrald& sp, std::vector<double>& coord)
{
	project(sp,coord);
	reproject_coordinates_to_positive(coord);
}


double spectrum_projector_3D::dot_product  (std::vector<double>& c1 , std::vector<double>& c2)
{
	assert(c1.size()==3 && c2.size()==3);
	return    m_A[0][0]*c1[0]*c2[0] + m_A[0][1]*c1[0]*c2[1] + m_A[0][2]*c1[0]*c2[2]
			+ m_A[1][0]*c1[1]*c2[0] + m_A[1][1]*c1[1]*c2[1] + m_A[1][2]*c1[1]*c2[2]
			+ m_A[2][0]*c1[2]*c2[0] + m_A[2][1]*c1[2]*c2[1] + m_A[2][2]*c1[2]*c2[2];
}

double spectrum_projector_3D::distance_spectrum_to_spectrum  (std::vector<double>& c1 , std::vector<double>& c2)
{
	std::vector<double> d;

	d.push_back(c1[0]-c2[0]);
	d.push_back(c1[1]-c2[1]);
	d.push_back(c1[2]-c2[2]);

	return std::sqrt(this->dot_product(d,d));
}


void spectrum_projector_3D::reproject (std::vector<double>& coord, const unsigned int x, const unsigned int y)
{
	const unsigned int z = 3 - x - y;

	double M_inverse [2][2];
	const double det = m_A[x][x] * m_A[y][y] - m_A[x][y] * m_A[y][x];
	M_inverse[0][0] = m_A[y][y] / det;
	M_inverse[1][1] = m_A[x][x] / det;
	M_inverse[0][1] = - m_A[x][y] / det;
	M_inverse[1][0] = - m_A[y][x] / det;

	double b [2];
	b[0] = m_A[x][0]*coord[0] + m_A[x][1]*coord[1] + m_A[x][2]*coord[2];
	b[1] = m_A[y][0]*coord[0] + m_A[y][1]*coord[1] + m_A[y][2]*coord[2];

	coord[x] = M_inverse[0][0] * b [0] + M_inverse[0][1] * b [1];
	coord[y] = M_inverse[1][0] * b [0] + M_inverse[1][1] * b [1];

	coord[z] = 0;
}

void spectrum_projector_3D::reproject (std::vector<double>& coord, const unsigned int x)
{
	const unsigned int y = (x+1)%3;
	const unsigned int z = (y+1)%3;
	coord[x] += coord[y] * m_A[x][y] / m_A[x][x];
	coord[x] += coord[z] * m_A[x][z] / m_A[x][x];
	coord[y] = 0;
	coord[z] = 0;
}

} // end namespace Fourier
} // end namespace ASTex
