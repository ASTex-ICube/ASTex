#include <stdlib.h>
#include "ASTex/Sparse/dictionaryProcessor.h"
#include "ASTex/easy_io.h"
#include "ASTex/utils.h"
#include "ASTex/rpn_utils.h"
#include "ASTex/Sparse/benchmarker.h"

int main(int argc, char **argv)
{	
	if(argc < 3)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " <in_texture> <out_directory>" << std::endl;
		return EXIT_FAILURE;
	}

	for(int i=1; i<argc; ++i)
	{
		ASTex::ImageRGBd im_in;
		ASTex::IO::loadu8_in_01(im_in, "/home/nlutz/img/" + std::string(argv[i]));

		std::string name_file = ASTex::IO::remove_path(argv[i]);
		std::string name_noext = ASTex::IO::remove_ext(name_file);

		SparseBenchmarker sb;
		sb.setInput(im_in);
		sb.setInputName(name_noext);
		sb.setRoot("/home/nlutz/img/sparse/");
		sb.setLearningDirectory("benchmarked");
		sb.setNbResults(3);
		sb.setLoadDictionaryMode(false);
		sb.setNumberIterationsLearning(10);
		sb.setNumberIterationsSynthesising(10);

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
		ASTex::DictionaryProcessor<ASTex::ImageRGBd> &dp = sb.dictionaryProcessor();
		dp.setInput(im_in);
		dp.setNbAtoms(32);
		dp.setSparsity(8);
		dp.setPatchSize(4, 4);
		dp.setPatchOffset(2, 2);
		for(unsigned j=4; j<10; ++j)
		{
			unsigned p = unsigned(std::pow(double(2), double(j)/2.0));
			dp.setPatchSize(p, p);
			sb.generate();
		}
	}
	return 0;
}
