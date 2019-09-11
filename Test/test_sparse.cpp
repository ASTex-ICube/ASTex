#include <stdlib.h>
#include "ASTex/Sparse/dictionaryProcessor.h"
#include "ASTex/easy_io.h"
#include "ASTex/utils.h"
#include "ASTex/rpn_utils.h"
int main(int argc, char **argv)
{	
	if(argc < 2)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << "<in_texture> <out_directory>" << std::endl;
		return EXIT_FAILURE;
	}

	ASTex::ImageRGBd im_in;
	ASTex::IO::loadu8_in_01(im_in, argv[1]);

	std::string name_file = ASTex::IO::remove_path(argv[1]);
	std::string name_noext = ASTex::IO::remove_ext(name_file);

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

	bool attemptToLoadDictionary = true;
	bool showDebugMessages = true;
	unsigned nbIterationsLearning = 512;
	unsigned nbIterationsSynthesis = 512;
	unsigned maxAtomsNumber = 512;
	unsigned sparsity = 4;
	itk::Size<2> patchSize;
	patchSize[0] = 4;
	patchSize[1] = 4;
	itk::Index<2> patchOffset;
	patchOffset[0] = 1;
	patchOffset[1] = 1;
	itk::Size<2> outputSize;
	outputSize[0] = 0;
	outputSize[1] = 0;

	std::string ioPath = std::string(argv[2]) + "/" + name_noext
				+ "_learn" + std::to_string(nbIterationsLearning)
				+ "_" + std::to_string(maxAtomsNumber) + "atoms"
				+ "_s" + std::to_string(sparsity)
				+ "_ox" + std::to_string(patchOffset[0])
				+ "_oy" + std::to_string(patchOffset[1])
				+ "_px" + std::to_string(patchSize[0])
				+ "_py" + std::to_string(patchSize[1]);
	ASTex::create_directory(ioPath);

	ASTex::DictionaryProcessor<ASTex::ImageRGBd> dp;

	if(!(attemptToLoadDictionary && dp.load(ioPath)))
	{
		dp.setInput(im_in);
		dp.setMaxAtomsNb(maxAtomsNumber);
		dp.setShowDebugMessages(showDebugMessages);
		dp.setSparsity(sparsity);
		dp.setPatchSize(patchSize[0], patchSize[1]);
		dp.setPatchOffset(patchOffset[0], patchOffset[1]);
		dp.dictionaryLearning(nbIterationsLearning);
		dp.save(ioPath);
	}
	dp.saveVizualisableAtoms(ioPath);

	dp.external_imageMatcher_path = "/home/nlutz/Git/colour-transfer/matchHistogram.sh";
	dp.synthesize(	outputSize[0]==0 ? im_in.width() : outputSize[0],
					outputSize[1]==0 ? im_in.height() : outputSize[1],
					nbIterationsSynthesis==0 ? dp.weights().size() * dp.sparsity() * 2 : nbIterationsSynthesis);

	return 0;
}
