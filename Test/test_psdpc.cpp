#include "ASTex/rpn_utils.h"
#include "ASTex/CSN/tpsd.h"
#include <ostream>
#include <ctime>

typedef struct
{
	Eigen::Vector2d vectors[2];
} CyclePair;

using CycleMapType = std::map<std::string, CyclePair>;

/**
 * @brief loadCycles load a cycles file.
 * Format of the file (don't copy the = symbols):
 * [<texture name>					//texture name with no extension or directory
 * <denominator> <denominator>		//two denomionators of the first cycle (cycle is [1/denominator]^2)
 * <denominator> <denominator>]*	//two denomionators of the second cycle
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
		if(!ifs.eof()) //I hate C++ sometimes
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


void psdToModulus(ImageSpectrald::PixelType &pix)
{
	pix = std::sqrt(pix);
}

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << "<in_texture> <out_texture>" << std::endl;
		return EXIT_FAILURE;
	}

	ImageRGBd input;
	std::vector<ImageGrayd> autocovariance;
	ImageGrayd output;
	ImageSpectrald ipsd;
	ImageSpectrald phase, randomPhase;

	std::string str_input = argv[1];
	std::string str_output = argv[2];
	std::string str_cycles = argv[3];
	create_directory(str_output);

	CyclePair cyclePair;

	std::string textureName = IO::remove_ext(IO::remove_path(str_input));
	if(argc < 4)
	{
		std::cerr << "Error: A cycle file must be provided to use the CSN with cycles!" << std::endl;
	}
	CycleMapType loadedCycles = loadCycles(str_cycles);
	CycleMapType::const_iterator cit = loadedCycles.find(textureName);
	if(cit == loadedCycles.end())
	{
		std::cerr << "Error: Texture name not found in the provided cycle map!" << std::endl;
		exit(EXIT_FAILURE);
	}
	cyclePair = (*cit).second;

	Eigen::Vector2d t0 = cyclePair.vectors[0], t1 = cyclePair.vectors[1];

	IO::loadu8_in_01(input, str_input);

	unsigned nbSamplesX = t0[0]*input.width(), nbSamplesY = t1[1]*input.height();

	//===============================//
	//DELETE THIS IDIOT

	int xOffset = 1024/21;
	int yOffset = 0;
	ImageRGBd idiot;
	idiot.initItk(std::round(t0[0]*input.width()), std::round(t1[1]*input.height()));
	for(unsigned nbTilesX=0; nbTilesX<unsigned(std::round(1.0/t0[0])); ++nbTilesX)
	{
		for(unsigned nbTilesY=0; nbTilesY<unsigned(std::round(1.0/t1[1])); ++nbTilesY)
		{
			idiot.for_all_pixels([&] (ImageRGBd::PixelType &pix, int x, int y)
			{
				pix = input.pixelAbsolute(unsigned((t0[0]*nbTilesX + t1[0]*nbTilesY)*input.width() + x + xOffset)%input.width(), unsigned((t1[1]*nbTilesY + t0[1]*nbTilesX)*input.height() + y + yOffset)%input.height());
			});
			IO::save01_in_u8(idiot, std::string("/home/nlutz/csn_tiles/") + textureName + "_" + std::to_string(nbTilesX) + "_" + std::to_string(nbTilesY) + ".png");
		}
	}


	//===============================//


//	autocovariance.resize(nbSamplesX*nbSamplesY);
//	ipsd.initItk(input.width(), input.height());
//	phase.initItk(input.width(), input.height());
//	randomPhase.initItk(input.width(), input.height());
//	output.initItk(input.width(), input.height());
//	Fourier::randomPhase(randomPhase);
//	for(unsigned x=0; x<nbSamplesX; ++x)
//	{
//		for(unsigned y=0; y<nbSamplesY; ++y)
//		{
//			ImageGrayd &ac = autocovariance[x*nbSamplesY+y];
//			Eigen::Vector2d t;
//			t[0] = x/double(nbSamplesX) * t0[0];
//			t[1] = y/double(nbSamplesY) * t1[1];
//			ac.initItk(1.0/t0[0], 1.0/t1[1], true);
//			ac.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x0, int x1)
//			{ //need to time average
//						pix += bilinear_interpolation(input, (t[0])*input.width(), (t[1])*input.height(), true)
//								* bilinear_interpolation(input, (t[0] + x0*t0[0] + x0*t1[0])*input.width(), (t[1]+x1*t1[1] + x1*t0[1])*input.height(), true);
//			});
//			Fourier::fftForwardModulusAndPhase(ac, ipsd, phase);
//			//IO::save01_in_u8(ac, str_output + "/ac_" + std::to_string(x) + "_" + std::to_string(y) + ".png");
//			//IO::save01_in_u8(ipsd, str_output + "/ipsd_" + std::to_string(x) + "_" + std::to_string(y) + ".png");
//			//ipsd.for_all_pixels(psdToModulus);
//			Fourier::fftInverseModulusAndPhase(ipsd, randomPhase, ac);
//			//IO::save01_in_u8(ac, str_output + "/rpn_" + std::to_string(x) + "_" + std::to_string(y) + ".png");
//			ac.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x0, int x1)
//			{
//				output.pixelAbsolute(x + (x0*t0[0]+x0*t1[0])*input.width(), y + (x1*t1[1]+x1*t0[1])*input.height()) = pix;
//			});
//		}
//	}
//	output.for_all_pixels([&] (ImageGrayd::PixelType &pix)
//	{
//			pix = pix > 1.0 ? 1.0 : (pix < 0.0 ? 0.0 : pix);
//	});
//	IO::save01_in_u8(output, "/home/nlutz/output.png");

	return 0;
}

