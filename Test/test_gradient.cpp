#include <cstdlib>

#include "ASTex/rpn_utils.h"
#include "ASTex/easy_io.h"

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cerr << "Usage: <scalar field image file> <output file> [16-bit?]" << std::endl;
		std::cerr << argv[0] << "" << std::endl;
		return EXIT_FAILURE;
	}

	std::string name_file = ASTex::IO::remove_path(argv[1]);
	std::string name_noext = ASTex::IO::remove_ext(name_file);

	using ImageTypeInput = ASTex::ImageGrayd;
	using ImageTypeOutput = ASTex::ImageRGBd;
	ImageTypeInput im_input;
	bool is16bit=false;
	if(argc>=4)
		is16bit = bool(atoi(argv[3]));

	if(!is16bit)
		IO::loadu8_in_01(im_input, argv[1]);
	else
	{
		ImageGrayu16 im_input16;
		im_input16.load(argv[1]);
		im_input.initItk(im_input16.width(), im_input16.height());
		im_input.for_all_pixels([&] (ImageTypeInput::PixelType &pix, int x, int y)
		{
			pix = im_input16.pixelAbsolute(x, y)/65535.0;
			//std::cout << pix << std::endl;
		});
	}
	ImageTypeOutput im_output = computeGradient(im_input);

	double minXpos=1.0, minYpos=1.0, maxXpos=0.5, maxYpos=0.5;
	double minXneg=0.5, minYneg=0.5, maxXneg=0, maxYneg=0;
	im_output.for_all_pixels([&] (ASTex::ImageRGBd::PixelType &pix)
	{
		if(pix[0]<0.5)
		{
			minXneg = std::min(minXneg, pix[0]);
			maxXneg = std::max(maxXneg, pix[0]);
		}
		if(pix[0]>0.5)
		{
			minXpos = std::min(minXpos, pix[0]);
			maxXpos = std::max(maxXpos, pix[0]);
		}
		if(pix[1]<0.5)
		{
			minYneg = std::min(minYneg, pix[1]);
			maxYneg = std::max(maxYneg, pix[1]);
		}
		if(pix[1]>0.5)
		{
			minYpos = std::min(minYpos, pix[1]);
			maxYpos = std::max(maxYpos, pix[1]);
		}
	});
	im_output.for_all_pixels([&] (ASTex::ImageRGBd::PixelType &pix)
	{
		if(pix[0]<0.5)
			pix[0] = ((pix[0]-minXneg)/(maxXneg - minXneg))/2.0;
		if(pix[0]>0.5)
			pix[0] = ((pix[0]-minXpos)/(maxXpos - minXpos))/2.0 + 0.5;
		if(pix[1]<0.5)
			pix[1] = ((pix[1]-minYneg)/(maxYneg - minYneg))/2.0;
		if(pix[1]>0.5)
			pix[1] = ((pix[1]-minYpos)/(maxYpos - minYpos))/2.0 + 0.5;
	});

	if(is16bit)
	{
		ASTex::ImageRGBu16 im_output16;
		im_output16.initItk(im_output.width(), im_output.height());
		im_output16.for_all_pixels([&] (ASTex::ImageRGBu16::PixelType &pix, int x, int y)
		{
			//std::cout << im_output.pixelAbsolute(x, y) << std::endl;
			for(unsigned i=0; i<3; ++i)
			{
				pix[i]=uint16_t(im_output.pixelAbsolute(x, y)[i]*65535.0);
			}
		});
		im_output16.save(argv[2]);
	}
	else
	{
		ASTex::IO::save01_in_u8(im_output, argv[2]);
	}

	return 0;
}
