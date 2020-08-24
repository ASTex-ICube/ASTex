#include "ASTex/rpn_utils.h"
#include "ASTex/CSN/tpsd.h"
#include <ostream>
#include <ctime>

void psdToModulus(ImageSpectrald::PixelType &pix)
{
	pix = std::sqrt(pix);
}

int trpn_mapmod(int argc, char **argv)
{
	if(argc < 4)
	{
		std::cerr << "Usage: " << argv[0]
				  << "<input psd> <input psd> <scalar map> <output texture..?>" << std::endl;
		exit(EXIT_FAILURE);
	}

	ImageSpectrald input, input2;
	ImageGrayd map;
	ImageGrayd output, output2;

	std::string str_input = argv[1];
	std::string str_input2 = argv[2];
	std::string str_map = argv[3];
	std::string str_output = argv[4];
	IO::loadu8_in_01(input, str_input);
	IO::loadu8_in_01(input2, str_input2);
	input.for_all_pixels(psdToModulus);
	input2.for_all_pixels(psdToModulus);
	IO::loadu8_in_01(map, str_map);

	ImageSpectrald phase;
	phase.initItk(input.width(), input.height());
	Fourier::randomPhase(phase);

	ImageSpectrald psd;
	psd.initItk(input.width(), input.height());

	Fourier::fftInverseModulusAndPhase(input, phase, output);
	Fourier::fftInverseModulusAndPhase(input2, phase, output2);
	double t0 = 0.2, t1 = 0.1;
	output.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
	{
		double dx = (double(x)/output.width())/t0;
		dx = dx - std::floor(dx); //fract
		dx = dx * (map.width()-1);

		double dy = (double(y)/output.height())/t1;
		dy = dy - std::floor(dy);
		dy = dy * (map.height()-1);

		ImageGrayd::PixelType w = bilinear_interpolation(map, dx, dy, true);
		psd.for_all_pixels([&] (ImageSpectrald::PixelType &pix, int x, int y)
		{
			pix = w*input.pixelAbsolute(x, y) + (1.0-w)*input2.pixelAbsolute(x, y);
		});

		Fourier::fftInverseModulusAndPhase(input, phase, output);

		pix = output.pixelAbsolute(x, y);
		std::cout << pix << std::endl;
	});

	output.for_all_pixels([] (ImageGrayd::PixelType &pix)
	{
		pix = std::min(std::max(0.0, 0.5 + 0.5), 1.0);
	});

	IO::save01_in_u8(output, str_output);

	return 0;
}

int trpn_cosmod(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cerr << "Usage: " << argv[0]
				  << "<input psd> <input psd> <output texture..?>" << std::endl;
		exit(EXIT_FAILURE);
	}

	ImageSpectrald input, input2;
	ImageGrayd output, output2;

	std::string str_input = argv[1];
	std::string str_input2 = argv[2];
	std::string str_output = argv[3];
	IO::loadu8_in_01(input, str_input);
	IO::loadu8_in_01(input2, str_input2);
	input.for_all_pixels(psdToModulus);
	input2.for_all_pixels(psdToModulus);

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

	output.for_all_pixels([] (ImageGrayd::PixelType &pix)
	{
		pix = std::min(std::max(0.0, 0.5 + 0.8*pix), 1.0);
	});

	IO::save01_in_u8(output, str_output);

	return 0;
}

int trpn_tpsd(int argc, char **argv)
{
	if(argc < 5)
	{
		std::cerr << "Usage: " << argv[0]
				  << "<input psd> <input psd> <input psd again> <input psd once more> <output texture..?>" << std::endl;
		exit(EXIT_FAILURE);
	}

	ImageSpectrald input, input2, input3, input4;
	ImageGrayd output, output2, output3, output4;

	std::string str_input = argv[1];
	std::string str_input2 = argv[2];
	std::string str_input3 = argv[3];
	std::string str_input4 = argv[4];
	std::string str_output = argv[5];

	IO::loadu8_in_01(input, str_input);
	IO::loadu8_in_01(input2, str_input2);
	IO::loadu8_in_01(input3, str_input3);
	IO::loadu8_in_01(input4, str_input4);

	input = scale(input, 800, 800);
	input2 = scale(input2, 800, 800);
	input3 = scale(input3, 800, 800);
	input4 = scale(input4, 800, 800);
	input.for_all_pixels(psdToModulus);
	input2.for_all_pixels(psdToModulus);
	input3.for_all_pixels(psdToModulus);
	input4.for_all_pixels(psdToModulus);

	unsigned nbSamplesX=320, nbSamplesY=160;

	CSN::TPSD<ImageSpectrald> tpsd;
	tpsd.setSize(input.width(), input.height());
	tpsd.setNbSamples(nbSamplesX, nbSamplesY);
	tpsd.setParallelogramCheat(true, false);
	std::function<std::function<void (ImageSpectrald::PixelType&, int, int)> (double, double)> func;
	func = [&] (double dx, double dy) -> std::function<void (ImageSpectrald::PixelType&, int, int)>
	{
		std::function<void (ImageSpectrald::PixelType&, int, int)> nestedFunc;
		nestedFunc = [&] (ImageSpectrald::PixelType &pix, int x, int y)
		{
			if(dx < 0.5)
				if(dy < 0.5)
				{
					double tdx = 2*dx, tdy = 2*dy;
					pix = input.pixelAbsolute(+x, y) * (1.0-tdx)*(1.0-tdy)
							+ input2.pixelAbsolute(x, y)*tdx*(1.0-tdy)
							+ input3.pixelAbsolute(x, y)*(1.0-tdx)*tdy
							+ input4.pixelAbsolute(x, y)*tdx*tdy;
				}
				else
				{
					double tdx = 2*dx, tdy = 2*(dy-0.5);
					pix = input3.pixelAbsolute(x, y) * (1.0-tdx)*(1.0-tdy)
							+ input4.pixelAbsolute(x, y)*tdx*(1.0-tdy)
							+ input.pixelAbsolute(x, y)*(1.0-tdx)*tdy
							+ input2.pixelAbsolute(x, y)*tdx*tdy;
				}
			else
			{
				if(dy < 0.5)
				{
					double tdx = 2*(dx-0.5), tdy = 2*dy;
					pix = input2.pixelAbsolute(x, y) * (1.0-tdx)*(1.0-tdy)
							+ input.pixelAbsolute(x, y)*tdx*(1.0-tdy)
							+ input4.pixelAbsolute(x, y)*(1.0-tdx)*tdy
							+ input3.pixelAbsolute(x, y)*tdx*tdy;
				}
				else
				{
					double tdx = 2*(dx-0.5), tdy = 2*(dy-0.5);
					pix = input4.pixelAbsolute(x, y) * (1.0-tdx)*(1.0-tdy)
							+ input3.pixelAbsolute(x, y)*tdx*(1.0-tdy)
							+ input2.pixelAbsolute(x, y)*(1.0-tdx)*tdy
							+ input.pixelAbsolute(x, y)*tdx*tdy;
				}
			}
		};
		return nestedFunc;
	};
	tpsd.setOperatorFunction(func);
	std::cout << "It passed the TPSD check" << std::endl;

	ImageSpectrald phase;
	phase.initItk(input.width(), input.height());
	Fourier::randomPhase(phase);

	output.initItk(input.width(), input.height());
	for(unsigned x=0; x<nbSamplesX; ++x)
		for(unsigned y=0; y<nbSamplesY; ++y)
		{
			ImageGrayd rpn;
			rpn.initItk(input.width(), input.height());
			ImageSpectrald psd = tpsd(x, y);
			psd.for_all_pixels(psdToModulus);
			Fourier::fftInverseModulusAndPhase(psd, phase, rpn);
			for(int yOut=y; yOut<output.height(); yOut+=nbSamplesY)
				for(int xOut=x; xOut<output.width(); xOut+=nbSamplesX)
					output.pixelAbsolute(xOut, yOut) = rpn.pixelAbsolute(xOut, yOut);
		}
	output.for_all_pixels([] (ImageGrayd::PixelType &pix)
	{
		pix = std::min(std::max(0.0, pix + 0.3), 1.0);
	});

	IO::save01_in_u8(output, str_output);
	return 0;
}

int main(int argc, char **argv)
{
	return trpn_tpsd(argc, argv);
}

