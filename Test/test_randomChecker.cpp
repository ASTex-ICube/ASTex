#include <stdlib.h>
#include "ASTex/easy_io.h"
#include <random>

using namespace ASTex;

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cerr << "Usage: <imageSize> <cellSize>" << std::endl;
		std::cerr << argv[0] << "" << std::endl;
		return EXIT_FAILURE;
	}

	unsigned checkerSize = atoi(argv[1]);
	unsigned cellSize = atoi(argv[2]);
	ImageRGBd checker;
	checker.initItk(checkerSize, checkerSize);

	std::default_random_engine generator;
	std::normal_distribution<double> distribution1(0.1, 0.01);
	std::normal_distribution<double> distribution2(0.9, 0.01);
	std::uniform_int_distribution<> distribInt(0, 1);

	for(unsigned chunkX = 0; chunkX<checkerSize/cellSize; ++chunkX)
	{
		for(unsigned chunkY = 0; chunkY<checkerSize/cellSize; ++chunkY)
		{
			ImageRGBd::PixelType pix;
			if(distribInt(generator)==0)
			{
				pix[0] = distribution1(generator);
				for(unsigned i=0; i<3; ++i)
				{
					pix[i] = pix[0];
				}
			}
			else
			{
				pix[0] = distribution2(generator);
				for(unsigned i=0; i<3; ++i)
				{
					pix[i] = pix[0];
				}
			}

			for(unsigned x=0; x<cellSize; ++x)
			{
				for(unsigned y=0; y<cellSize; ++y)
				{
					checker.pixelAbsolute(chunkX*cellSize+x, chunkY*cellSize+y) = pix;
				}
			}
		}
	}

	checker.for_all_pixels([&] (ImageRGBd::PixelType &pix)
	{
		for(unsigned i=0; i<3; ++i)
		{
			pix[i] = std::min(std::max(0.0, pix[i]), 1.0);
		}
	});

	IO::save01_in_u8(checker, ${ASTEX_TEMPO_PATH}+"checkerG.png");

	return 0;
}
