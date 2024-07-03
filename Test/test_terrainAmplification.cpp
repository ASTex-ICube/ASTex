#include <stdlib.h>
#include <iostream>
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
		std::cerr << argv[0] << " <im_source> <im_target> [arg list]" << std::endl;
		return EXIT_FAILURE;
	}

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

	std::string name_file = ASTex::IO::remove_path(argv[1]);
	std::string name_noext = ASTex::IO::remove_ext(name_file);

	ASTex::ImageRGBd im_source, im_target, im_output;
	IO::loadu8_in_01(im_source, argv[1]);
	IO::loadu8_in_01(im_target, argv[2]);
	ASTex::DictionaryProcessor<ASTex::ImageRGBd> *dp = new ASTex::DictionaryProcessor<ASTex::ImageRGBd>;

	unsigned numberIterationsLearning = 30;
	unsigned nbAtoms = 180;
	unsigned sparsity = 2;
	unsigned patchOffset = 4;
	unsigned patchSize = 12;
	unsigned amplificationFactor = 3; //should ideally be

	if(argc >= 4)
	{
		std::ifstream inputArgs(argv[3]);
		if(inputArgs)
		{
			inputArgs >> numberIterationsLearning >> nbAtoms >> sparsity >> patchOffset >> patchSize >> amplificationFactor;
		}
		inputArgs.close();
	}

	std::string ioPath = std::string(TEMPO_PATH+"dictionaries") + "/" + name_noext
			+ "_learn" + std::to_string(numberIterationsLearning)
			+ "_" + std::to_string(nbAtoms) + "atoms"
			+ "_s" + std::to_string(sparsity)
			+ "_ox" + std::to_string(patchOffset)
			+ "_p" + std::to_string(patchSize);
	if(!dp->load(ioPath))
	{
		dp->setInput(im_source);
		dp->setMaxAtomsNb(nbAtoms);
		dp->setPatchOffset(patchOffset, patchOffset);
		dp->setPatchSize(patchSize, patchSize);
		dp->setSparsity(sparsity);
		dp->setShowDebugMessages(true);
		dp->dictionaryLearning(numberIterationsLearning);
		dp->save(ioPath);
		dp->saveVizualisableAtoms(ioPath);
	}


	ASTex::ImageRGBd test_output = dp->reconstructInput();
	IO::save01_in_u8(test_output, TEMPO_PATH+"terrain7Reconstructed.png");

	im_output = dp->amplify<CompareRGBSum>(im_target, amplificationFactor);
	IO::save01_in_u8(im_output, TEMPO_PATH+"terrainAmplified.png");

	return 0;
}
