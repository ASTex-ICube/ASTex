#include "ASTex/rpn_utils.h"
#include <ostream>
#include <ctime>

//Generates the spectrum of an image

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cerr << "Usage: " << argv[0]
				  << "<input psd> <input psd> <output texture..?>" << std::endl;
		exit(EXIT_FAILURE);
	}

	ImageSpectrald input;
	ImageSpectrald input2;
	ImageGrayd output;
	ImageGrayd output2;

	std::string str_input = argv[1];
	std::string str_input2 = argv[2];
	std::string str_output = argv[3];

	IO::loadu8_in_01(input, str_input);
	IO::loadu8_in_01(input2, str_input2);

	output.initItk(input.width(), input.height());
	output2.initItk(input.width(), input.height());
	ImageSpectrald phase;
	phase.initItk(input.width(), input.height());
	Fourier::randomPhase(phase);

	Fourier::fftInverseModulusAndPhase(input, phase, output);
	Fourier::fftInverseModulusAndPhase(input2, phase, output2);
	double f0 = 7.0, f1 = 3.0;
	double phi = M_PI;
	double energy1=1.0, energy2=2.0;
	output.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
	{
		pix *= 0.5 * (	energy1 * std::sqrt(2) * cos(2*M_PI*f0*(double(x)/output.width())) * cos(2*M_PI*f1*y/output.height())
					+	energy2 * std::sqrt(2) * cos(2*M_PI*f0*(double(x)/output.width()) + phi) * cos(2*M_PI*f1*y/output.height() + phi)) + 0.5;
		pix = std::min(std::max(0.0, pix), 1.0);
		std::cout << pix << std::endl;
	});

	IO::save01_in_u8(output, str_output);
}
