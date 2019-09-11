#include "ASTex/rpn_utils.h"
#include <ostream>

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cerr << "Usage: " << argv[0]
				  << "<output directory> <texture file 1> [texture file 2] ... [texture file n]" << std::endl;
		exit(EXIT_FAILURE);
	}

	ImageRGBd input;
	ImageRGBd output;

	std::string out_dir = argv[1];

	for(unsigned i=2; i<argc; ++i)
	{
		std::string input_path = argv[i];
		std::string name_file = IO::remove_path(argv[i]);
		std::string name_noext = IO::remove_ext(name_file);
		create_directory(out_dir);

		unsigned seed = getpid();
		std::cout << "Using pid as seed: " << seed << std::endl;
		srand(seed);

		IO::loadu8_in_01(input, input_path);
		colored_RPN(input, output);

		output.for_all_pixels([&] (ImageRGBd::PixelType &pix)
		{
			for(unsigned i=0; i<3; ++i)
			{
				pix[i] = std::max(std::min(1.0, pix[i]), 0.0);
			}
		});
		IO::save01_in_u8(output, out_dir + "/" + name_noext + "_RPN_" + std::to_string(seed) + ".png");
	}
}
