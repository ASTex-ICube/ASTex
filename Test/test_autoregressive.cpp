#include <stdlib.h>
#include <ASTex/autoregressive.h>
#include <ASTex/easy_io.h>

int main(int argc, char **argv)
{
	if(argc < 1)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " <in_texture>" << std::endl;
		return EXIT_FAILURE;
	}

	Autoregressive<ImageGrayd> ar;
	double mean = 0;
	double variance = 0.1;
	ar.setWhiteNoiseParameters(mean, variance);
	ar.setOrder(3, 3);
	ar.setCoeff(-3, 0, -0.5);
	ar.setCoeff(-2, 0, 0.5);
	ar.setCoeff(-1, 0, -0.5);

	ar.setConstant(0.5);

	ar.setWidth(64);
	ar.setHeight(64);
	ImageGrayd im_out = ar.simulateFromLeftUp();
	im_out.for_all_pixels([&] (ImageGrayd::PixelType &pix)
	{
		pix = std::max(std::min(1.0, pix), 0.0);
	});
	IO::save01_in_u8(im_out, "/home/nlutz/im_ar.png");
	return 0;
}
