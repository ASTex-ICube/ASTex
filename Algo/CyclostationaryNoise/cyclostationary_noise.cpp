#include <stdlib.h>
#include "ASTex/easy_io.h"
#include "ASTex/CSN/csn_texture.h"
#include "ASTex/image_rgb.h"
#include "ASTex/Stamping/stamper.h"
#include <map>
#include "ASTex/histogram.h"

using namespace ASTex;

typedef struct
{
	Eigen::Vector2d vectors[2];
} CyclePair;

using CycleMapType = std::map<std::string, CyclePair>;

/**
 * @brief loadCycles load a cycles file.
 * Format of the file (don't copy the = symbols):
 * [<texture name>					//texture name with no extension or directory
 * <denominator> <denominator>		//two denomionators of the first cycle
 * <denominator> <denominator>]*	//two denomionators of the second cycle
 * if denominator is an integer, then it is interpreted as 1/denominator (useful for periodic textures).
 * if denominator is a floating point between 0 and 1, then it is interpreted as is.
 * =======
 * Exemple:
 * =======
 * bricks5
 * 14	20
 * 7	0
 *
 * herringbone2
 * 7	0
 * 0	26
 * =======
 */
CycleMapType loadCycles(std::string filename)
{
	CycleMapType knownCycles;
	std::string name;
	CyclePair cycles;
	std::ifstream ifs(filename);
	assert(ifs);
	while(!ifs.eof())
	{
		ifs >> name;
		if(!ifs.eof())
		{
			double readNumber1, readNumber2;
            ifs >> readNumber1 >> readNumber2;
			if((readNumber1>0 && readNumber1<1) || (readNumber2>0 && readNumber2<1))
				cycles.vectors[0] = Eigen::Vector2d(readNumber1, readNumber2);
			else
				cycles.vectors[0] = Eigen::Vector2d(readNumber1 == 0 ? 0 : 1.0/readNumber1,
													readNumber2 == 0 ? 0 : 1.0/readNumber2);
			ifs >> readNumber1 >> readNumber2;
			if((readNumber1>0 && readNumber1<1) || (readNumber2>0 && readNumber2<1))
				cycles.vectors[1] = Eigen::Vector2d(readNumber1, readNumber2);
			else
				cycles.vectors[1] = Eigen::Vector2d(readNumber1 == 0 ? 0 : 1.0/readNumber1,
													readNumber2 == 0 ? 0 : 1.0/readNumber2);
			knownCycles.insert({name, cycles});
		}
	}
	return knownCycles;
}

typedef struct
{
	bool useCycles;
	double gamma;
	unsigned outputWidth, outputHeight;
	unsigned proceduralBlendingMode;
	bool usePca;
	bool useGaussianTransfer;
	bool useYCbCr;
	bool useCyclicTransfer;
	double uvScale;
	double cyclicTransferRadius;
	unsigned cyclicTransferSamples;
} ArgumentsType;

/**
 * @brief loadArguments load an argument file.
 * Format of the file (don't copy the = symbols):
 * =======
 * <use cycles?>
 * <gamma>							//weight's exponant (for HPN)
 * <output width> <output height>
 * <procedural blending mode>		//determines the blending mode (0: HPN, 1: spot noise)
 * <use PCA>
 * <use Gaussian transfer>
 * <use YCbCr>
 * <use cyclic transfer (needs gaussian transfer and cycles)>
 * <uvScale>
 * <transfer radius>
 * <transfer number of samples>
 * =======
 * Example:
 * =======
 * 1
 * 2.0
 * 4096 4096
 * 0
 * 0
 * 1
 * 1
 * 1
 * 0.8
 * 1.5
 * 9
 * =======
 */
ArgumentsType loadArguments(std::string filename)
{
	ArgumentsType arguments;
	std::ifstream ifs(filename);
	assert(ifs);
	unsigned uValue;
	ifs >> uValue;
	arguments.useCycles = uValue;
	ifs >> arguments.gamma;
	ifs >> arguments.outputWidth >> arguments.outputHeight;
	ifs >> arguments.proceduralBlendingMode;
	ifs >> uValue;
	arguments.usePca = uValue;
	ifs >> uValue;
	arguments.useGaussianTransfer = uValue;
	ifs >> uValue;
	arguments.useYCbCr = uValue;
	ifs >> uValue;
	arguments.useCyclicTransfer = uValue;
	ifs >> arguments.uvScale;
	ifs >> arguments.cyclicTransferRadius;
	ifs >> arguments.cyclicTransferSamples;
	return arguments;
}

int main(int argc, char **argv)
{
	std::setlocale(LC_NUMERIC, "fr_FR");
	if(argc < 3)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " <in_texture> <out_texture> <argument file> [cycle file]" << std::endl;
		return EXIT_FAILURE;
	}

	using ImageType = ImageRGBd;
	using PcaImageType = CSN::CSN_Texture<ImageType>::PcaImageType;

	ImageType im_in;
	std::string filename_in = std::string(argv[1]);
	std::string out_filename = std::string(argv[2]);
	std::string filename_arguments = std::string(argv[3]);
	std::string filename_cycles = std::string(argv[4]);
	IO::loadu8_in_01(im_in, filename_in);
	std::string textureName = IO::remove_ext(IO::remove_path(filename_in));

	ArgumentsType arguments = loadArguments(filename_arguments);

	CyclePair cyclePair;
	if(arguments.useCycles)
	{
		if(argc < 4)
		{
			std::cerr << "Error: A cycle file must be provided to use the CSN with cycles!" << std::endl;
		}
		CycleMapType loadedCycles = loadCycles(filename_cycles);
		CycleMapType::const_iterator cit = loadedCycles.find(textureName);
		if(cit == loadedCycles.end())
		{
			std::cerr << "Error: Texture name not found in the provided cycle map!" << std::endl;
			exit(EXIT_FAILURE);
		}
		cyclePair = (*cit).second;
	}

	auto textonProcedure = [&] (const PcaImageType &image) -> PcaImageType
	{
		PcaImageType centeredImage;
		centeredImage.initItk(image.width(), image.height());
		centeredImage.copy_pixels(image);

		HistogramRGBBase<typename PcaImageType::DataType> histo(centeredImage);
		PcaImageType::PixelType mean = histo.meanPixelType();
		PcaImageType meanImage;
		meanImage.initItk(image.width(), image.height(), true);
		meanImage.for_all_pixels([&] (typename PcaImageType::PixelType &pix, int x, int y)
		{
			if(arguments.useCycles)
			{
				unsigned cxIndex = 0;
				unsigned cyIndex = 1; //need to find this automatically because it's not always true.
				unsigned maxCyclesX = unsigned(std::round(1.0/cyclePair.vectors[cxIndex][0]));
				unsigned maxCyclesY = unsigned(std::round(1.0/cyclePair.vectors[cyIndex][1]));
				for(unsigned cx=0; cx<maxCyclesX; ++cx)
				{
					for(unsigned cy=0; cy<maxCyclesY; ++cy)
					{
						double xShift= x	+ (cx*cyclePair.vectors[cxIndex][0]
											+ cy*cyclePair.vectors[cyIndex][0])*image.width();
						double yShift= y	+ (cx*cyclePair.vectors[cxIndex][1]
											+ cy*cyclePair.vectors[cyIndex][1])*image.height();
						pix += bilinear_interpolation(image, xShift, yShift, true);
					}
				}
				pix = pix * (1.0/(maxCyclesX*maxCyclesY));
			}
			else
				pix = mean;
		});

		//Histogram<ImageRGBd>::loadImageFromCsv(meanImage, "/home/nlutz/mean.csv");
		centeredImage.for_all_pixels([&] (PcaImageType::PixelType &pix, int x, int y)
		{
			pix[0] -= meanImage.pixelAbsolute(x, y)[0];
			pix[1] -= meanImage.pixelAbsolute(x, y)[1];
			pix[2] -= meanImage.pixelAbsolute(x, y)[2];
			pix = pix * (1.0/std::sqrt(centeredImage.width()*centeredImage.height()));
		});
		double alpha = 0.02;
		//normalization
		std::cout << "applying smooth transition function" << std::endl;
		PcaImageType::PixelType norm = PcaImageType::zero(), normPhi = PcaImageType::zero();
		centeredImage.for_all_pixels([&] (PcaImageType::PixelType &pix, int x, int y)
		{
			for(unsigned i=0; i<3; ++i)
			{
				norm[i] += pix[i] * pix[i];
			}
		});
		for(unsigned i=0; i<3; ++i)
		{
			norm[i] = sqrt(norm[i]);
		}
//		ImageGrayd d;
//		d.initItk(centeredImage.width(), centeredImage.height());
//		d.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
//		{
//			double distanceToBorder;
//			double distanceX = double(std::min(x, centeredImage.width()-1-x));
//			double distanceY = double(std::min(y, centeredImage.height()-1-y));
//			distanceToBorder = std::min(distanceX/centeredImage.width(), distanceY/centeredImage.height());
//			if(distanceToBorder<alpha)
//			{
//				for(unsigned i=0; i<3; ++i)
//				{
//					pix = sqrt(distanceToBorder/alpha);
//				}
//			}
//			else
//				pix = 1.0;
////			std::cout << pix << std::endl;
//		});
		//IO::save01_in_u8(d, "/home/nlutz/d.png");
		centeredImage.for_all_pixels([&] (PcaImageType::PixelType &pix, int x, int y)
		{
			double distanceToBorder;
			double distanceX = double(std::min(x, centeredImage.width()-1-x));
			double distanceY = double(std::min(y, centeredImage.height()-1-y));
			distanceToBorder = std::min(distanceX/centeredImage.width(), distanceY/centeredImage.height());
			for(unsigned i=0; i<3; ++i)
			{
				pix[i] = pix[i]/norm[i];
				if(distanceToBorder<alpha)
				{
					pix[i] = pix[i] * sqrt(distanceToBorder/alpha);
				}
				normPhi[i] += pix[i]*pix[i];
			}
		});
		for(unsigned i=0; i<3; ++i)
		{
			normPhi[i] = sqrt(normPhi[i]);
		}
		centeredImage.for_all_pixels([&] (PcaImageType::PixelType &pix, int x, int y)
		{
			for(unsigned i=0; i<3; ++i)
			{
				pix[i] /= normPhi[i];
				pix[i] *= norm[i];
			}
		});

		Stamping::StampDiscrete<PcaImageType> stamp(centeredImage);
		stamp.setInterpolationRule(Stamping::StampDiscrete<PcaImageType>::BILINEAR_PERIODIC);
		if(arguments.useCycles)
		{
			srand(0);
			Stamping::SamplerCycles sampler;
			sampler.setNbPoints(420);
			sampler.setCycles(cyclePair.vectors[0].cast<float>()/2.0, cyclePair.vectors[1].cast<float>()/2.0);
			Stamping::StamperTexton<PcaImageType> stamper(&sampler, &stamp);
			stamper.setPeriodicity(true);
			stamper.setUseMargins(false);
			stamper.setSpot(false);
			centeredImage = stamper.generate(	arguments.outputWidth == 0 ? im_in.width()*2 : arguments.outputWidth,
												arguments.outputHeight == 0 ? im_in.height()*2 : arguments.outputHeight);
		}
		else
		{
			Stamping::SamplerUniform sampler;
			sampler.setNbPoints(1000);
			Stamping::StamperTexton<PcaImageType> stamper(&sampler, &stamp);
			stamper.setPeriodicity(true);
			stamper.setSpot(false);
			centeredImage = stamper.generate(	arguments.outputWidth == 0 ? im_in.width()*1: arguments.outputWidth,
												arguments.outputHeight == 0 ? im_in.height()*1 : arguments.outputHeight);
		}
		centeredImage.for_all_pixels([&] (PcaImageType::PixelType &pix, int x, int y)
		{
			itk::Index<2> index;
			index[0] = x%meanImage.width();
			index[1] = y%meanImage.height();
			for(int i=0; i<3; ++i)
			{
				pix[i] += meanImage.pixelAbsolute(index)[i];
			}
		});
		IO::save01_in_u8(meanImage, "/home/nlutz/average.png");
		return centeredImage;
	};

	CSN::CSN_Texture<ImageType> csn;
	csn.setTexture(im_in);
	//ImageRGBd cycleEvaluationMap = csn.debug_cycleEvaluationMap(129, 129, Eigen::Vector2d(cyclePair.vectors[0][0], cyclePair.vectors[1][1]), 0.1);
	//IO::save01_in_u8(cycleEvaluationMap, std::string("/home/nlutz/cycleEvaluationMap129_") + textureName + ".png");
	csn.setCycles(cyclePair.vectors[0], cyclePair.vectors[1]);
//	std::cout << csn.testCycles_v2(im_in, cyclePair.vectors[0], cyclePair.vectors[1], 40, 40) << std::endl;
	csn.setUseCycles(arguments.useCycles);
	csn.setGamma(arguments.gamma);
	csn.setUsePca(arguments.usePca);
	csn.setUseGaussianTransfer(arguments.useGaussianTransfer);
	csn.setUseYCbCr(arguments.useYCbCr);
	csn.setUseCyclicTransfer(arguments.useCyclicTransfer);
	csn.setUVScale(arguments.uvScale);
	if(arguments.useCyclicTransfer)
		csn.setCyclicTransferPolicy(arguments.cyclicTransferRadius, arguments.cyclicTransferSamples);
	if(arguments.proceduralBlendingMode == 1)
	{
		csn.setProceduralBlendingSubstitute(textonProcedure);
	}
//	if(arguments.useCycles)
//	{
//		csn.estimateCycles(cyclePair.vectors[0], cyclePair.vectors[1], 0.02, true, 32);
//	}
	std::cout << textureName << std::endl;
	std::cout << "Proposed cycle x: " << std::endl << cyclePair.vectors[0] << std::endl;
	std::cout << "Proposed cycle y: " << std::endl << cyclePair.vectors[1] << std::endl;
	std::cout << "Estimated cycle x: " << std::endl << csn.cycleX() << std::endl;
	std::cout << "Estimated cycle y: " << std::endl << csn.cycleY() << std::endl;
	ImageType output = csn.synthesize(arguments.outputWidth, arguments.outputHeight);
	output.for_all_pixels([&] (ImageType::PixelType &pix)
	{
		for(int i=0; i<3; ++i)
		{
			pix[i] = pix[i] > 1.0 ? 1.0 : (pix[i] < 0.0 ? 0.0 : pix[i]);
		}
	});

	IO::save01_in_u8(im_in, "/home/nlutz/input.png");
	IO::save01_in_u8(output, out_filename);
	return 0;
}
