#include "ASTex/rpn_utils.h"
#include "ASTex/CSN/tpsd.h"
#include <ostream>
#include <ctime>

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

	ImageGrayd input;
	std::vector<ImageGrayd> autocovariance;
	ImageGrayd output;
	ImageSpectrald ipsd;
	ImageSpectrald phase, randomPhase;

	std::string str_input = argv[1];
	std::string str_output = argv[2];
	create_directory(str_output);

	Eigen::Vector2d t0({0.25, 0}), t1({0, 0.25});

	IO::loadu8_in_01(input, str_input);

	unsigned nbSamplesX = 4, nbSamplesY = 4;

	autocovariance.resize(nbSamplesX*nbSamplesY);
	ipsd.initItk(input.width(), input.height());
	phase.initItk(input.width(), input.height());
	randomPhase.initItk(input.width(), input.height());
	Fourier::randomPhase(randomPhase);
	for(unsigned x=0; x<nbSamplesX; ++x)
	{
		for(unsigned y=0; y<nbSamplesY; ++y)
		{
			ImageGrayd &ac = autocovariance[x*nbSamplesY+y];
			Eigen::Vector2d t;
			t[0] = x/double(nbSamplesX) * t0[0];
			t[1] = y/double(nbSamplesY) * t1[1];
			ac.initItk(input.width(), input.height(), true);
			ac.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x0, int x1)
			{ //need to time average
				for(unsigned k=0; k<4; ++k)
					for(unsigned l=0; l<4; ++l)
					{
						pix += bilinear_interpolation(input, (t[0]+k*t0[0])*input.width(), (t[1]+l*t1[1])*input.height(), true)
								* bilinear_interpolation(input, (t[0]+k*t0[0])*input.width() + x0, (t[1]+l*t1[1])*input.height() + x1, true);
					}
				pix = pix * (1.0/(4*4));
			});
			Fourier::fftForwardModulusAndPhase(ac, ipsd, phase);
			IO::save01_in_u8(ac, str_output + "/ac_" + std::to_string(x) + "_" + std::to_string(y) + ".png");
			IO::save01_in_u8(ipsd, str_output + "/ipsd_" + std::to_string(x) + "_" + std::to_string(y) + ".png");
			//ipsd.for_all_pixels(psdToModulus);
			Fourier::fftInverseModulusAndPhase(ipsd, randomPhase, ac);
			IO::save01_in_u8(ac, str_output + "/rpn_" + std::to_string(x) + "_" + std::to_string(y) + ".png");
		}
	}

	return 0;
}

