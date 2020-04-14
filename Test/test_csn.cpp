#include <stdlib.h>
#include "ASTex/easy_io.h"
#include "ASTex/CSN/csn_texture.h"
#include "ASTex/image_rgb.h"
#include <map>

using namespace ASTex;

typedef struct
{
	Eigen::Vector2d vectors[2];
} CyclePair;

using CycleMapType = std::map<std::string, CyclePair>;

/**
 * @brief loadCycles load a cycles file.
 * Format of the file (don't copy the = symbols):
 * [<texture name>					//texture name with no extension or directory
 * <denominator> <denominator>		//two denomionators of the first cycle (cycle is [1/denominator]^2)
 * <denominator> <denominator>]*	//two denomionators of the second cycle
 * =======
 * Exemple:
 * =======
 * bricks5
 * 14	20
 * 7	0
 *
 * herringbone2
 * 7	0
 * 0	26
 * =======
 */
CycleMapType loadCycles(std::string filename)
{
	CycleMapType knownCycles;
	std::string name;
	CyclePair cycles;
	std::ifstream ifs(filename);
	assert(ifs);
	while(!ifs.eof())
	{
		ifs >> name;
		if(!ifs.eof()) //I hate C++ sometimes
		{
            unsigned readNumber1, readNumber2;
            ifs >> readNumber1 >> readNumber2;
			cycles.vectors[0] = Eigen::Vector2d(readNumber1 == 0 ? 0 : 1.0/readNumber1,
												readNumber2 == 0 ? 0 : 1.0/readNumber2);
			ifs >> readNumber1 >> readNumber2;
			cycles.vectors[1] = Eigen::Vector2d(readNumber1 == 0 ? 0 : 1.0/readNumber1,
												readNumber2 == 0 ? 0 : 1.0/readNumber2);
			knownCycles.insert({name, cycles});
		}
	}
	return knownCycles;
}

typedef struct
{
	bool useCycles;
	double gamma;
	unsigned outputWidth, outputHeight;
} ArgumentsType;

/**
 * @brief loadArguments load an argument file.
 * Format of the file (don't copy the = symbols):
 * =======
 * <use cycles?>
 * <gamma>
 * <output width> <output height>
 * =======
 * Exemple:
 * =======
 * 1
 * 2.0
 * 4096 4096
 * =======
 */
ArgumentsType loadArguments(std::string filename)
{
	ArgumentsType arguments;
	std::ifstream ifs(filename);
	assert(ifs);
	unsigned uValue;
	ifs >> uValue;
	arguments.useCycles = uValue;
	ifs >> arguments.gamma;
	ifs >> arguments.outputWidth >> arguments.outputHeight;
	return arguments;
}

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " <in_texture> <out_texture> <argument file> [cycle file]" << std::endl;
		return EXIT_FAILURE;
	}

	using ImageType = ImageRGBd;

	ImageType im_in;
	std::string filename_in = std::string(argv[1]);
	std::string out_filename = std::string(argv[2]);
	std::string filename_arguments = std::string(argv[3]);
	std::string filename_cycles = std::string(argv[4]);
	IO::loadu8_in_01(im_in, filename_in);
	std::string textureName = IO::remove_ext(IO::remove_path(filename_in));

	ArgumentsType arguments = loadArguments(filename_arguments);

	CyclePair cyclePair;
	if(arguments.useCycles)
	{
		if(argc < 4)
		{
			std::cerr << "Error: A cycle file must be provided to use the CSN with cycles!" << std::endl;
		}
		CycleMapType loadedCycles = loadCycles(filename_cycles);
		CycleMapType::const_iterator cit = loadedCycles.find(textureName);
		if(cit == loadedCycles.end())
		{
			std::cerr << "Error: Texture name not found in the provided cycle map!" << std::endl;
			exit(EXIT_FAILURE);
		}
		cyclePair = (*cit).second;
	}

	CSN::CSN_Texture<ImageType> csn;
	csn.setTexture(im_in);
	csn.setCycles(cyclePair.vectors[0], cyclePair.vectors[1]);
	csn.setUseCycles(arguments.useCycles);
	csn.setGamma(arguments.gamma);
	ImageType output = csn.synthesize(arguments.outputWidth, arguments.outputHeight);
	IO::save01_in_u8(output, out_filename);
	return 0;
}
