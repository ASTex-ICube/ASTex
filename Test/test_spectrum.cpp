#include "ASTex/rpn_utils.h"
#include <ostream>
#include <ctime>

//Generates the spectrum of an image

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cerr << "Usage: " << argv[0]
				  << "<input texture> <output spectrum>" << std::endl;
		exit(EXIT_FAILURE);
	}

	ImageRGBd input;
	ImageRGBd output;

	std::string str_input = argv[1];
	std::string str_output = argv[2];

	IO::loadu8_in_01(input, str_input);

	ImageGrayd ch0, ch1, ch2;
	extract3Channels(input, ch0, ch1, ch2);

	ImageSpectrald sp0, sp1, sp2, phase;
	sp0.initItk(input.width(), input.height());
	sp1.initItk(input.width(), input.height());
	sp2.initItk(input.width(), input.height());
	phase.initItk(input.width(), input.height());

	Fourier::fftForwardModulusAndPhase(ch0, sp0, phase);
	Fourier::fftForwardModulusAndPhase(ch1, sp1, phase);
	Fourier::fftForwardModulusAndPhase(ch2, sp2, phase);

	fold3Channels(output, sp0, sp1, sp2);
	IO::save01_in_u8(output, str_output);
}
