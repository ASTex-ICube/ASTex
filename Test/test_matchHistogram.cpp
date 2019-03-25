#include <stdlib.h>
#include "ASTex/Sparse/dictionaryProcessor.h"
#include "ASTex/easy_io.h"
#include "ASTex/utils.h"
#include "ASTex/rpn_utils.h"

int main(int argc, char **argv)
{
	if(argc < 4)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " <in_texture> <in_model> <out_filename>" << std::endl;
		return EXIT_FAILURE;
	}

	using ImageType = ASTex::ImageRGBd;

	ImageType im_in;
	ASTex::IO::loadu8_in_01(im_in, std::string(argv[1]));

	ImageType im_model;
	ASTex::IO::loadu8_in_01(im_model, std::string(argv[2]));

	std::string out_filename = argv[3];

	matchImage(im_in, im_model);

	IO::save01_in_u8(im_in, out_filename);

	return 0;
}
