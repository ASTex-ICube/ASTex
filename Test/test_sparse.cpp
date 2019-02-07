#include <stdlib.h>
#include "ASTex/Sparse/dictionaryProcessor.h"
#include "ASTex/easy_io.h"
#include "ASTex/utils.h"
#include "ASTex/rpn_utils.h"

int main(int argc, char **argv)
{	
	if(argc < 3)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " <in_texture> <out_directory>" << std::endl;
		return EXIT_FAILURE;
	}

	ASTex::ImageRGBd im_in;
	ASTex::IO::loadu8_in_01(im_in, std::string(argv[1]));

	std::string out_dir = argv[2];
	std::string name_file = ASTex::IO::remove_path(argv[1]);
	std::string name_noext = ASTex::IO::remove_ext(name_file);
	ASTex::create_directory(out_dir);
	out_dir = out_dir + "/" + name_noext;
	ASTex::create_directory(out_dir);

	class CompareRGBSum
	{
	public:
		CompareRGBSum() {}
		bool operator()(const ASTex::ImageRGBd::PixelType &pix, const ASTex::ImageRGBd::PixelType &other)
		{
			double normpix = pix[0] + pix[1] + pix[2];
			double normother = other[0] + other[1] + other[2];
			return normpix < normother;
		}
	};
	ASTex::DictionaryProcessor<ASTex::ImageRGBd> dp;
	dp.setInput(im_in);
	dp.setNbAtoms(32);
	dp.setSparsity(8);
	dp.setPatchSize(6, 6);
	dp.setPatchOffset(2, 2);
	dp.readInput();
	dp.dictionaryLearning<CompareRGBSum>(4);

	for(unsigned i=0; i<dp.nbAtoms(); ++i)
	{
		ASTex::IO::save01_in_u8(dp.dictionary().atom(i).content(), out_dir + "/d" + std::to_string(i) + ".png");
	}

	ASTex::ImageRGBd reconstructedInput = dp.reconstructInput();
	ASTex::IO::save01_in_u8(reconstructedInput, "/home/nlutz/reconstructedInput.png");


	return 0;
}
