#include <stdlib.h>
#include "ASTex/easy_io.h"
#include "ASTex/CSN/csn_texture.h"
#include "ASTex/image_rgb.h"

using namespace ASTex;

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " <in_texture> <out_texture>" << std::endl;
		return EXIT_FAILURE;
	}

	ImageRGBd im_in;
	IO::loadu8_in_01(im_in, std::string(argv[1]));
	std::string out_filename = argv[2];

	CSN::CSN_Texture<ImageRGBd> csn;
	csn.setTexture(im_in);
	//csn.setCycles(Eigen::Vector2d(1.0/14.0, 1.0/20.0), Eigen::Vector2d(1.0/7.0, 0))
	csn.setCycles(Eigen::Vector2d(1.0/4.0, 0), Eigen::Vector2d(0, 1.0/4.0));
	csn.setUseCycles(true);
	ImageRGBd output = csn.synthesize(4096, 4096);
	IO::save01_in_u8(output, out_filename);
	return 0;
}
