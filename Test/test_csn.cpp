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
	ImageRGBd output = csn.synthesize(4096, 4096);
	IO::save01_in_u8(output, out_filename);
	return 0;
}
