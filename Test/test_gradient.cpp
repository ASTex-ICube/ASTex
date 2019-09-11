#include <stdlib.h>
#include "ASTex/easy_io.h"
#include "ASTex/rpn_utils.h"

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cerr << "Usage: <scalar field image file> <output directory>" << std::endl;
		std::cerr << argv[0] << "" << std::endl;
		return EXIT_FAILURE;
	}

	std::string name_file = ASTex::IO::remove_path(argv[1]);
	std::string name_noext = ASTex::IO::remove_ext(name_file);

	ASTex::ImageGrayd im_source;
	ASTex::IO::loadu8_in_01(im_source, argv[1]);

	ASTex::ImageRGBd im_output = computeGradient(im_source);
	ASTex::IO::save01_in_u8(im_output, argv[2]);

	return 0;
}
