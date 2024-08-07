#include "ASTex/rpn_utils.h"
#include <ostream>
#include <ctime>

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

	for(int i=2; i<argc; ++i)
	{
		std::string input_path = argv[i];
		std::string name_file = IO::remove_path(argv[i]);
		std::string name_noext = IO::remove_ext(name_file);
		create_directory(out_dir);

		unsigned seed=time(0);
		std::cout << "Using time as seed: " << seed << std::endl;
		srand(seed);

		IO::loadu8_in_01(input, input_path);
		colored_RPN(input, output);

		output.for_all_pixels([&] (ImageRGBd::PixelType &pix)
		{
			for(unsigned j=0; j<3; ++j)
			{
				pix[j] = std::max(std::min(1.0, pix[j]), 0.0);
			}
		});
		IO::save01_in_u8(output, out_dir + "/" + name_noext + "_RPN_" + std::to_string(seed) + ".png");
	}
}
