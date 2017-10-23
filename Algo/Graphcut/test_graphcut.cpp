#include "graphcut.h"

using namespace ASTex;

int main(int argc, char** argv)
{
	time_t t = time(NULL);
	//t = 10235234:
	//t = 1456765179;

	std::srand(t);
	std::cout << "Random number generator seed: " << t << std::endl;
	
	// Jelly2 with this seed:
	//std::srand(1456746461);
	// blocked 58 -> ...

	ImageRGBu8 inputImage;
	std::string inputName;

	if( argc < 2 )
		inputName=TEMPO_PATH+"single-leaf-photo-03.jpg";
	else
		inputName = std::string(argv[1]);

	// Load input image

	inputImage.load(inputName);

	int input_w = inputImage.width();
	int input_h = inputImage.height();
	
	std::cout << "Input image resolution: " << input_w << "*" << input_h << std::endl;

	size_t sep = inputName.find_last_of("\\/");
	if (sep != std::string::npos)
		inputName = inputName.substr(sep + 1, inputName.size() - sep - 1);

	std::string textureName;
	std::string extensionName;

	size_t dot = inputName.find_last_of(".");
	if (dot != std::string::npos)
	{
		textureName = inputName.substr(0, dot);
		extensionName = inputName.substr(dot, inputName.size() - dot);
	}
	else
	{
		textureName = inputName;
		extensionName = "";
	}

	std::cout << "Texture name: " << textureName << std::endl;

	std::string outputNameInitLoc = textureName + "_init_gcsynthesis.png";
	std::string outputNameInitLocSeams = textureName + "_init_s_gcsynthesis.png";

	std::string outputNameLoc = textureName + "_gcsynthesis.png";
	std::string outputNameLocSeams = textureName + "_s_gcsynthesis.png";
	std::cout << "Synthesis result image will be written as: " << outputNameLoc << std::endl;

	GCTexture gctex;

	gctex.initialize(inputImage, textureName, 1200, 1200);
	gctex.setSeamSize(1);
	gctex.setBlendingParams(5, 2);
	gctex.setSSDStep(3);
	gctex.topDownRandomFilling(input_w / 3, input_h / 3, false);

	std::cout << "Initial synthesis done!\n";
	gctex.writeOutputImage(outputNameInitLoc, outputNameInitLocSeams);

	gctex.writeOutputImage(outputNameLoc, outputNameLocSeams);
	std::cout << "Output image written\n";

	return EXIT_SUCCESS;
}

