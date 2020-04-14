#include <stdlib.h>
#include "ASTex/easy_io.h"
#include "ASTex/rpn_utils.h"
#include "ASTex/mipmap.h"

int main(int argc, char **argv)
{
	if(argc < 5)
	{
		std::cerr << "Usage: <output filename> <input heightmap> <tile reduction factor> <num. of tiles> [<texture for layer i> <max altitude (0-1) of layer i>]+" << std::endl;
		std::cerr << argv[0] << "" << std::endl;
		return EXIT_FAILURE;
	}

	std::string output_filename = argv[1];
//	std::string name_texture = ASTex::IO::remove_path(argv[2]);
//	std::string name_noext = ASTex::IO::remove_ext(name_texture);

	ASTex::ImageGrayd heightmap;
	ASTex::IO::loadu8_in_01(heightmap, argv[2]);

	ASTex::ImageRGBAd colorMap;

	unsigned tileReductionFactor = std::atoi(argv[3]);
	unsigned numberOfTiles = std::atoi(argv[4]);

	std::vector<std::pair<Mipmap<ImageRGBd>, double>> textureVector; //something else than a vector is better but w/e

	for(int i=5; i<argc; i+=2)
	{
		ASTex::ImageRGBd texture;
		ASTex::IO::loadu8_in_01(texture, argv[i]);
		Mipmap<ASTex::ImageRGBd> mipmap(texture);
		mipmap.setMode(ISOTROPIC);
		mipmap.setMaxPowReductionLevel(tileReductionFactor);
		mipmap.generate();

		double maxAltitude = std::atof(argv[i+1]);
		textureVector.push_back(std::make_pair(mipmap, maxAltitude));
	}


	//integrity check
	int width = textureVector[0].first.texture().width();
	int height = textureVector[0].first.texture().height();
	double maxAltitude = textureVector[0].second;
	for(auto &texturePair : textureVector)
	{
		if(width != texturePair.first.texture().width())
		{
			std::cerr << "Error: not all textures have the same width" << std::endl;
			exit(EXIT_FAILURE);
		}
		if(height != texturePair.first.texture().height())
		{
			std::cerr << "Error: not all textures have the same height" << std::endl;
			exit(EXIT_FAILURE);
		}
		if(maxAltitude>texturePair.second || texturePair.second>1 || texturePair.second<0)
		{
			std::cerr << "Error: altitudes must be all be sorted from lowest to highest and between 0 and 1" << std::endl;
			exit(EXIT_FAILURE);
		}
		maxAltitude = texturePair.second;
	}
	textureVector.back().second = 1.0;

	//let's get to work

	double widthScale = double(heightmap.width()-1)/width;
	double heightScale = double(heightmap.height()-1)/height;

	colorMap.initItk(width*numberOfTiles, height*numberOfTiles);
	colorMap.for_all_pixels([&] (ASTex::ImageRGBAd::PixelType &pix, int x, int y)
	{
		unsigned i=0;
		for(; bilinear_interpolation(heightmap, (x*widthScale)/numberOfTiles, heightmap.height() - (y*heightScale)/numberOfTiles - 1, false)>textureVector[i].second; ++i);
		ASTex::ImageRGBd &texture = textureVector[i].first.mipmap(tileReductionFactor, tileReductionFactor);
		pix[0] = texture.pixelAbsolute(x%texture.width(), y%texture.height())[0];
		pix[1] = texture.pixelAbsolute(x%texture.width(), y%texture.height())[1];
		pix[2] = texture.pixelAbsolute(x%texture.width(), y%texture.height())[2];
		pix[3] = 0.75;
	});

	ASTex::IO::save01_in_u8(colorMap, output_filename);

	return 0;
}
