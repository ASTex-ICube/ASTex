#include "graphcut.h"

using namespace ASTex;


//
//

int main(int argc, char** argv)
{
	time_t t = time(NULL);

	std::srand(t);
	std::cout << "Random number generator seed: " << t << std::endl;
	
	// Jelly2 with this seed:
	//std::srand(1456746461);
	// blocked 58 -> ...

	ImageRGBu8 inputImage;

	if (argc < 2)
	{
		std::cerr << "Usage: test_graphcut input_image.png\n";
		return EXIT_FAILURE;
	}

	std::string inputName(argv[1]);

	// Load input image

	inputImage.load(inputName);

//	typedef itk::ImageFileReader< ImageRGBu8::ItkImg > ReaderLocalType;
//	ReaderLocalType::Pointer reader = ReaderLocalType::New();
//	reader->SetFileName(inputName);
//	try
//	{
//		reader->Update();
//	}
//	catch (itk::ExceptionObject & error)
//	{
//		std::cerr << "Failed to load file " << inputName << std::endl;
//		std::cerr << error << std::endl;
//		return EXIT_FAILURE;
//	}
//	inputImage.itk() = reader->GetOutput();

	int input_w = inputImage.width();
	int input_h = inputImage.height();
	
	std::cout << "Input image resolution: " << input_w << "*" << input_h << std::endl;

	size_t sep = inputName.find_last_of("\\/");
	if (sep != std::string::npos)
		inputName = inputName.substr(sep + 1, inputName.size() - sep - 1);

	std::string textureName;
//	std::string extensionName;

	size_t dot = inputName.find_last_of(".");
	if (dot != std::string::npos)
	{
		textureName = inputName.substr(0, dot);
//		extensionName = inputName.substr(dot, inputName.size() - dot);
	}
	else
	{
		textureName = inputName;
//		extensionName = "";
	}

	std::cout << "Texture name: " << textureName << std::endl;

	std::string outputNameInitLoc = textureName + "_init_gcsynthesis.png";
	std::string outputNameInitLocSeams = textureName + "_init_s_gcsynthesis.png";

	std::string outputNameLoc = textureName + "_gcsynthesis.png";
	std::string outputNameLocSeams = textureName + "_s_gcsynthesis.png";
	std::cout << "Synthesis result image will be written as: " << outputNameLoc << std::endl;

	GCTexture gctex;

	gctex.initialize(inputImage, textureName, 275, 275); // 275
	gctex.setSeamSize(1);
	gctex.setBlendingParams(5, 2);
	gctex.setSSDStep(1);
	gctex.topDownRandomFilling(input_w / 3, input_h / 3, false);
	std::cout << "Initial synthesis done!\n";
	gctex.writeOutputImage(outputNameInitLoc, outputNameInitLocSeams);

	gctex.setBorderSize(5 + 3);
	for (int k = 0; k < 100; k++)
	{
		gctex.seamsErrorSubPatchRefinement(100, 3, false);
	}

	// Blending enabled
	std::cout << "Refinement done!\n";
	gctex.writeOutputImage(outputNameLoc, outputNameLocSeams);
	std::cout << "Output image written\n";

	return EXIT_SUCCESS;
}

