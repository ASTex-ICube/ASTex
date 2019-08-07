#ifndef __NOISE__H__
#define __NOISE__H__

#include <cmath>
#include <ASTex/special_io.h>
#include <ASTex/fourier.h>
#include <ASTex/local_spectrum.h>
#include <ASTex/utils.h>
#include <ASTex/colorspace_filters.h>
#include <ASTex/mask.h>
#include <ASTex/distances_maps.h>
#include <ASTex/pca.h>
#include <ASTex/easy_io.h>
#include "histogram.h"

using namespace ASTex;

template <int N>
inline double log_scale(double d)
{
	const double x = 1.0 / N;
	return (std::log(d + x) - std::log(x)) / (std::log(1 + x) - std::log(x));
}

class CompareWeightedPCA
{
public:
	CompareWeightedPCA() {}
	bool operator()(const typename ImageRGBd::PixelType &object, const typename ImageRGBd::PixelType &other)
	{
		return object[0] * weight[0] + object[1] * weight[1] + object[2] * weight[2]
			<  other[0]  * weight[0] + other[1]  * weight[1] + other[2]  * weight[2];
	}
	static double weight[3];
};

typedef enum{ALL=0, NORMAL=1, ONE_PHASE=2, SAME_PHASE=3, NULL_PHASE=4, SEVERAL_PHASES=5, TEXTON=6, PCA_SYNTHESIS=7, LUMINANCE=8} color_dephasing_mode_t;
static const std::string color_dephasing_mode_map[]={"all", "normal", "one_phase", "same_phase", "null_phase", "several_phases", "texton", "pca", "luminance"};

typedef enum{RGB_SPACE=0, LAB_SPACE=1, LUV_SPACE=2} color_space_t;
static const std::string color_space_map[]={"RGB", "Lab", "Luv"};

void rpn_function(const std::string &out_path, const std::string& in_texture, color_dephasing_mode_t mode=NORMAL, int doubledScale=1);

template <typename REAL, typename std::enable_if<std::is_floating_point<REAL>::value>::type* = nullptr >
void saveHistogram2csv(CommonSpectral<ImageGrayBase<REAL>, false>& spectrum, const std::string &filename)
{
	std::set<REAL> histogram;

	spectrum.for_all_pixels([&histogram] (REAL& pix)
	{
	   histogram.insert(pix);
	});

	std::ofstream csvFile;

	try
	{
		csvFile.open (filename);
		std::for_each(histogram.begin(), histogram.end(), [&csvFile] (REAL value)
		{
			csvFile << value << std::endl;
		});
		csvFile.close();
	}
	catch(const std::exception& ex)
	{
		std::cerr << "Error while writting csv file : " << ex.what() << std::endl;
	}

	return;
}

template <typename REAL, typename std::enable_if<std::is_floating_point<REAL>::value>::type* = nullptr >
REAL volumeSpectral(const CommonSpectral<ImageGrayBase<REAL>, false>& modulus)
{
	REAL volume=REAL(0);

	modulus.for_all_pixels([&volume] (REAL& pix)
	{
		volume+=pix;
	});

	return volume/sqrt(modulus.width()*modulus.height());
}

void saveFourierModulusPhaseGray(const std::string &out_path, const std::string& in_texture);

void rpn_scalar(const ImageSpectrald& modulus, ImageSpectrald& phase, ImageGrayd& output);

template <typename REAL, typename std::enable_if<std::is_floating_point<REAL>::value>::type* = nullptr >
void extractXYZ(const ImageCommon<ImageRGBBase<REAL>, false>& rgbImage,
				ImageCommon<ImageRGBBase<REAL>, false>& xyzImage)
{
	xyzImage.initItk(rgbImage.width(), rgbImage.height());

	xyzImage.parallel_for_all_pixels([&rgbImage] (typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType& pix, int x, int y)
	{
		pix=ColorSpace::colorRGBtoXYZ<  typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType,
										typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType>
			(rgbImage.pixelAbsolute(x, y));
	});

	return;
}

template <typename REAL, typename std::enable_if<std::is_floating_point<REAL>::value>::type* = nullptr >
void extractLuv(const ImageCommon<ImageRGBBase<REAL>, false>& rgbImage,
				ImageCommon<ImageGrayBase<REAL>, false>& LChannel,
				ImageCommon<ImageGrayBase<REAL>, false>& uChannel,
				ImageCommon<ImageGrayBase<REAL>, false>& vChannel)
{
	LChannel.initItk(rgbImage.width(), rgbImage.height());
	uChannel.initItk(rgbImage.width(), rgbImage.height());
	vChannel.initItk(rgbImage.width(), rgbImage.height());

	rgbImage.parallel_for_all_pixels([&LChannel, &uChannel, &vChannel] (const typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType& pix, int x, int y)
	{
		typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType LuvPix
				= ColorSpace::colorXYZtoLUV<  typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType,
											  typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType>
					(ColorSpace::colorRGBtoXYZ<  typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType,
												 typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType>
					 (pix));

		LChannel.pixelAbsolute(x, y)=LuvPix[0];
		uChannel.pixelAbsolute(x, y)=LuvPix[1];
		vChannel.pixelAbsolute(x, y)=LuvPix[2];
	});

	return;
}

template <typename REAL, typename std::enable_if<std::is_floating_point<REAL>::value>::type* = nullptr >
void extractLab(const ImageCommon<ImageRGBBase<REAL>, false>& rgbImage,
				ImageCommon<ImageGrayBase<REAL>, false>& LChannel,
				ImageCommon<ImageGrayBase<REAL>, false>& aChannel,
				ImageCommon<ImageGrayBase<REAL>, false>& bChannel)
{
	LChannel.initItk(rgbImage.width(), rgbImage.height());
	aChannel.initItk(rgbImage.width(), rgbImage.height());
	bChannel.initItk(rgbImage.width(), rgbImage.height());

	rgbImage.parallel_for_all_pixels([&LChannel, &aChannel, &bChannel] (const typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType& pix, int x, int y)
	{
		typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType LabPix
				= ColorSpace::colorXYZtoLAB<  typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType,
											  typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType>
					(ColorSpace::colorRGBtoXYZ<  typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType,
												 typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType>
					 (pix));

		LChannel.pixelAbsolute(x, y)=LabPix[0];
		aChannel.pixelAbsolute(x, y)=LabPix[1];
		bChannel.pixelAbsolute(x, y)=LabPix[2];
	});

	return;
}

template <typename REAL, typename std::enable_if<std::is_floating_point<REAL>::value>::type* = nullptr >
void extract3Channels(const ImageCommon<ImageRGBBase<REAL>, false>& rgbImage,
				ImageCommon<ImageGrayBase<REAL>, false>& redChannel,
				ImageCommon<ImageGrayBase<REAL>, false>& greenChannel,
				ImageCommon<ImageGrayBase<REAL>, false>& blueChannel)
{
	redChannel.initItk(rgbImage.width(), rgbImage.height());
	greenChannel.initItk(rgbImage.width(), rgbImage.height());
	blueChannel.initItk(rgbImage.width(), rgbImage.height());

	unsigned int channel;

	auto channel2gray = [&channel, &rgbImage] (typename ImageCommon<ImageGrayBase<REAL>, false>::PixelType& pix, int x, int y)
	{
		pix=rgbImage.pixelAbsolute(x, y)[channel];
	};

    channel=0;
	redChannel.for_all_pixels(channel2gray);

    channel=1;
	greenChannel.for_all_pixels(channel2gray);

    channel=2;
	blueChannel.for_all_pixels(channel2gray);
	return;
}

template <typename REAL, typename std::enable_if<std::is_floating_point<REAL>::value>::type* = nullptr >
void foldRGBfromXYZ(ImageCommon<ImageRGBBase<REAL>, false>& rgbImage,
			 const ImageCommon<ImageRGBBase<REAL>, false>& xyzImage)
{
	if(!rgbImage.is_initialized() || rgbImage.width()!=xyzImage.width() || rgbImage.height()!=xyzImage.height())
		rgbImage.initItk(xyzImage.width(), xyzImage.height());

	rgbImage.parallel_for_all_pixels([&xyzImage] (typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType& pix, int x, int y)
	{
		pix=ColorSpace::colorXYZtoRGB<  typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType,
										typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType>
			(xyzImage.pixelAbsolute(x, y));
	});
}

template <typename REAL, typename std::enable_if<std::is_floating_point<REAL>::value>::type* = nullptr >
void foldRGBfromLuv(ImageCommon<ImageRGBBase<REAL>, false>& rgbImage,
			 const ImageCommon<ImageRGBBase<REAL>, false>& LuvImage)
{
	if(!rgbImage.is_initialized() || rgbImage.width()!=LuvImage.width() || rgbImage.height()!=LuvImage.height())
		rgbImage.initItk(LuvImage.width(), LuvImage.height());

	rgbImage.parallel_for_all_pixels([&LuvImage] (typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType& pix, int x, int y)
	{
		pix=ColorSpace::colorXYZtoRGB<  typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType,
										typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType>
			(ColorSpace::colorLUVtoXYZ< typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType,
										typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType>
			(LuvImage.pixelAbsolute(x, y)));
	});
}

template <typename REAL, typename std::enable_if<std::is_floating_point<REAL>::value>::type* = nullptr >
void foldRGBfromLab(ImageCommon<ImageRGBBase<REAL>, false>& rgbImage,
			 const ImageCommon<ImageRGBBase<REAL>, false>& LabImage)
{
	if(!rgbImage.is_initialized() || rgbImage.width()!=LabImage.width() || rgbImage.height()!=LabImage.height())
		rgbImage.initItk(LabImage.width(), LabImage.height());

	rgbImage.parallel_for_all_pixels([&LabImage] (typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType& pix, int x, int y)
	{
		pix=ColorSpace::colorXYZtoRGB<  typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType,
										typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType>
			(ColorSpace::colorLABtoXYZ< typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType,
										typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType>
			(LabImage.pixelAbsolute(x, y)));
	});
}

template <typename REAL, typename std::enable_if<std::is_floating_point<REAL>::value>::type* = nullptr >
void fold3Channels(ImageCommon<ImageRGBBase<REAL>, false>& channeledImage,
			 const ImageCommon<ImageGrayBase<REAL>, false>& channel1,
			 const ImageCommon<ImageGrayBase<REAL>, false>& channel2,
			 const ImageCommon<ImageGrayBase<REAL>, false>& channel3)
{
	int width = channel1.width();
	int height = channel1.height();
	assert(channel2.width() == width && channel3.width() == width);
	assert(channel2.height() == height && channel3.height() == height);
    if(!channeledImage.is_initialized() || channeledImage.width()!=channel1.width() || channeledImage.height()!=channel1.height())
        channeledImage.initItk(channel1.width(), channel1.height());

	channeledImage.for_all_pixels([&channel1, &channel2, &channel3] (typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType& pix, int x, int y)
    {
        pix[0]=channel1.pixelAbsolute(x, y);
        pix[1]=channel2.pixelAbsolute(x, y);
		pix[2]=channel3.pixelAbsolute(x, y);
    });
}

template <typename REAL, typename std::enable_if<std::is_floating_point<REAL>::value>::type* = nullptr >
void image2periodicComponent(const ImageCommon<ImageGrayBase<REAL>, false>& input,
							 ImageCommon<ImageGrayBase<REAL>, false>& periodicComponent)
{
	if(!input.is_initialized())
		return;

	int W = input.width();
	int H = input.height();

	if(!periodicComponent.is_initialized() || periodicComponent.width()!=W || periodicComponent.height()!=H)
		periodicComponent.initItk(W, H);

	// 1 - compute boundary intensity (v)
	ImageGrayd bound(W,H,true);
	for (int x=0;x<W;x++)
	{
		double v = input.pixelAbsolute(x,0) - input.pixelAbsolute(x,H-1);
		bound.pixelAbsolute(x,0) += v;
		bound.pixelAbsolute(x,H-1) -= v;
	}
	for (int y=0;y<H;y++)
	{
		double v = input.pixelAbsolute(0,y) - input.pixelAbsolute(W-1,y);
		bound.pixelAbsolute(0,y) += v;
		bound.pixelAbsolute(W-1,y) -= v;
	}

	// 2 - DFT of v
	using  FFTType = itk::ForwardFFTImageFilter< ImageGrayd::ItkImg, ImageSpectralcd::ItkImg >;
	FFTType::Pointer fftFilter = FFTType::New();
	fftFilter->SetInput(bound.itk());
	fftFilter->Update();
	ImageSpectralcd ft_bound(fftFilter->GetOutput());

	// 3 - divide by cosin....
	//
	// compute coef by row
	double cx = 2.0*M_PI/double(W);
	std::vector<double> coef_x(W);
	for (int x=0; x<W; ++x)
		coef_x[x] = 2*std::cos(cx*x);
	// by column
	double cy = 2.0*M_PI/double(H);
	std::vector<double> coef_y(H);
	for (int y=0; y<H; ++y)
		coef_y[y] = 2*std::cos(cy*y);
	// apply computations on all pixels
	ft_bound.for_all_pixels([&] (ImageSpectralcd::PixelType& P, int x, int y)
	{
		P /= 4 - coef_x[x] - coef_y[y]; // WARNING OPPOSITE SIGN COMPARE TO ARTICLE
	});
	// except for (0,0) -> 0,0
	ft_bound.pixelAbsolute(0,0) = ImageSpectralcd::PixelType(0,0);

	// 4 - FFT inverse => smooth
	//
	using  IFFTType = itk::InverseFFTImageFilter< ImageSpectralcd::ItkImg, ImageGrayd::ItkImg >;
	IFFTType::Pointer ifftFilter = IFFTType::New();
	ifftFilter->SetInput(ft_bound.itk());
	ifftFilter->Update();
	ImageGrayd smooth(ifftFilter->GetOutput());

	// 4' - compute periodic = input - smooth
	//
	for(auto its = std::make_tuple(periodicComponent.beginIterator(),input.beginConstIterator(),smooth.beginConstIterator());
		 !std::get<0>(its).IsAtEnd(); ++std::get<0>(its), ++std::get<1>(its), ++std::get<2>(its))
	{
		std::get<0>(its).Value() = std::get<1>(its).Value() - std::get<2>(its).Value();
	}

	return;
}

template <typename REAL, typename std::enable_if<std::is_floating_point<REAL>::value>::type* = nullptr >
void image2periodicComponent(const ImageCommon<ImageRGBBase<REAL>, false>& input,
							 ImageCommon<ImageRGBBase<REAL>, false>& periodicComponent)
{
	if(!input.is_initialized())
		return;

	if(!periodicComponent.is_initialized() || periodicComponent.width()!=input.width() || periodicComponent.height()!=input.height())
		periodicComponent.initItk(input.width(), input.height());

	ImageCommon<ImageGrayBase<REAL>, false> redChannel, greenChannel, blueChannel;
	extract3Channels(input, redChannel, greenChannel, blueChannel);

	image2periodicComponent(redChannel, redChannel);
	image2periodicComponent(greenChannel, greenChannel);
	image2periodicComponent(blueChannel, blueChannel);

	fold3Channels(periodicComponent, redChannel, greenChannel, blueChannel);
	return;
}

template <typename REAL, typename std::enable_if<std::is_floating_point<REAL>::value>::type* = nullptr >
void fourierModulus_color(const CommonSpectral<ImageGrayBase<REAL>, false>& modulus, ImageCommon<ImageRGBBase<REAL>, false> &output)
{
	if(!modulus.is_initialized())
		return;

	if(!output.is_initialized() || output.width()!=modulus.width() || output.height()!=modulus.height())
		output.initItk(modulus.width(), modulus.height());

	REAL min0 = compute_min(modulus);
	REAL max0 = compute_max(modulus);

	output.for_all_pixels([&] (typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType &pix, int x, int y)
	{
		typename ImageGrayBase<REAL>::PixelType q = modulus.pixelAbsolute(x, y);
		double q2 = (q - min0) / (max0 - min0);
		double p = log_scale<100000>(q2);

		//construct hsv pixel
		typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType hsvPixel;
		hsvPixel[0]=120-p*120;
		hsvPixel[1]=255;
		hsvPixel[2]=255;

		pix = ColorSpace::colorHSVtoRGB<typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType,
										typename ImageCommon<ImageRGBBase<REAL>, false>::PixelType>(hsvPixel);
	});

	return;
}

void gray_RPN(const ImageGrayd& in, ImageGrayd& out, unsigned int extendX, unsigned int extendY,  bool crop=false, bool periodic_component=true, bool call_srand=true);

void colored_RPN(const ImageRGBd& in, ImageRGBd& out, color_space_t colorSpace=RGB_SPACE, color_dephasing_mode_t mode=NORMAL, unsigned int extendX=0, unsigned int extendY=0,
				 bool crop=false, bool periodic_component=true, bool call_srand=false, double scale_randomPhase=1.0,
				 const ImageSpectrald *phase=NULL);

double compute_crossCorrelation_diff(const ImageRGBd& in1, const ImageRGBd& in2, int channel1, int channel2, ImageGrayd& diff);

void matchImage(ImageRGBd& source, const ImageRGBd& target, std::string external_program_name = std::string());

template<typename NUMTYPE>
void matchImage(ImageCommon<ImageGrayBase<NUMTYPE>, false>& image, const ImageCommon<ImageGrayBase<NUMTYPE>, false>& input)
{
	using ImageType = ImageCommon<ImageGrayBase<NUMTYPE>, false>;
	using HistogramType = HistogramGrayBase<NUMTYPE>;
	HistogramType histInput(input), cdfImage(image), histHistoryOfImage(image);
	unsigned count = 0;
	for(typename HistogramType::iterator itThis=cdfImage.begin(); itThis!=cdfImage.end(); ++itThis)
	{ //transforms the histogram into a cdf
		count += (*itThis).second;
		(*itThis).second = count;
	}

	//compute translations from a sorted version of image to image
	ImageRGB32 colorMap;
	colorMap.initItk(input.width(), input.height());
	colorMap.for_all_pixels([&] (ImageRGB32::PixelType &pix, int x, int y)
	{
		int position = (*cdfImage.find(image.pixelAbsolute(x, y))).second
				- (*histHistoryOfImage.find(image.pixelAbsolute(x, y))).second--;
		pix[0] = position%colorMap.width(), pix[1] = position/colorMap.width();
	});


	ImageType orderedInput;
	orderedInput.initItk(input.width(), input.height());
	typename HistogramType::iterator itInput = histInput.begin();
	orderedInput.for_all_pixels([&] (typename ImageType::PixelType &pix)
	{
		if((*itInput).second == 0)
			++itInput;
		pix = (*itInput).first;
		--(*itInput).second;
	});

	ImageType matchedImage;
	matchedImage.initItk(image.width(), image.height());
	matchedImage.for_all_pixels([&] (typename ImageType::PixelType &pix, int x, int y)
	{
//		std::cout << "(" << x << ", " << y << ") -> " << "x: " << colorMap.pixelAbsolute(x, y)[0] << ", y: " << colorMap.pixelAbsolute(x, y)[1] << std::endl;
		matchedImage.pixelAbsolute(	x, y ) = orderedInput.pixelAbsolute(colorMap.pixelAbsolute(x, y)[0], colorMap.pixelAbsolute(x, y)[1]);
	});

	image = matchedImage;
}

ImageRGBd computeGradient(const ImageGrayd &heightField);

#endif
