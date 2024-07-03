#include <stdlib.h>
#include <ctime>
#include "ASTex/ContentExchange/multiTextureProcessor.h"

int main(int argc, char **argv)
{
	if(argc < 4)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " <out_directory> <in_tile> <argument file> [[in_texture i] [periodicity of texture i]]*" << std::endl;
		return EXIT_FAILURE;
	}

	//argument file.csv:
	//nb patches
	//nb contents per patch
	//outputWidth outputHeight
	//seed (0 to call srand(time(NULL)))
	//#Transformation list, as a bunch of 2x2 matrices, columns first.
	//#for example, rotate by 180 degrees and scale of 1.1 will be:
	//-1.1 0 0 -1.1
	//#If you want to have either a rotation or a scale or both, you have to be exhaustive in your argument list.

	using ImageType = ImageRGBd;

	ImageRGBd im_tile;
	IO::loadu8_in_01(im_tile, std::string(argv[2]));

	std::string out_dir = argv[1];
	std::string name_file = IO::remove_path(argv[2]);
	std::string name_noext = IO::remove_ext(name_file);
	create_directory(out_dir);

	ContentExchange::MultiTextureProcessor<ImageType> pProcessor(im_tile);
	ASTex::TexturePool<ImageType, double> texturePool;

	///Argument file
	std::ifstream ifs_args(argv[3]);
	if(!ifs_args)
	{
		std::cerr << "Error: could not open argument file" << std::endl;
		exit(EXIT_FAILURE);
	}
	unsigned nb_patches, nb_contents, output_width, output_height, seed, gradientManipulation;
	ifs_args >> nb_patches >> nb_contents >> output_width >> output_height >> seed >> gradientManipulation;
	while(!ifs_args.eof())
	{
		double m00, m01, m02, m10, m11, m12, m20, m21, m22;
		ifs_args >> m00;
		if(!ifs_args.eof())
		{
			ifs_args >> m01 >> m02 >> m10 >> m11 >> m12 >> m20 >> m21 >> m22;
			Eigen::Matrix3d M;
			M(0, 0)=m00;
			M(0, 1)=m01;
			M(0, 2)=m02;
			M(1, 0)=m10;
			M(1, 1)=m11;
			M(1, 2)=m12;
			M(2, 0)=m20;
			M(2, 1)=m21;
			M(2, 2)=m22;
			texturePool.addTransformationMatrix(M);
		}
	}
	texturePool.addTexture(im_tile, true);
	texturePool.setTexturesAreGradient(gradientManipulation);

	for(int i=4; i<argc - (argc%2==0 ? 0 : 2); i+=2)
	{
		ImageType im_texture;
		IO::loadu8_in_01(im_texture, argv[i]);
		unsigned periodicity = atoi(argv[i+1]);
		texturePool.addTexture(im_texture, periodicity);
	}
	if(argc%2!=0)
	{
		ImageType im_texture;
		IO::loadu8_in_01(im_texture, argv[argc-1]);
		texturePool.addTexture(im_texture, false);
	}
	texturePool.generate();

	pProcessor.setSeed(seed == 0 ? time(nullptr) : seed);

	pProcessor.setFilteringMode(ISOTROPIC);
	pProcessor.setNbContentsPerPatch(nb_contents);
	pProcessor.setOutputSize(output_width, output_height);
	pProcessor.patches_initRandom(nb_patches);
	pProcessor.patches_dilate(false);
	pProcessor.contents_initDefault();
	pProcessor.setTexturePool(texturePool);
	pProcessor.contents_explorePools();

//	for(unsigned p=0; p<pProcessor.nbPatches(); ++p)
//		for(unsigned c=0; c<pProcessor.nbContents(); ++c)
//			IO::save01_in_u8(pProcessor.patchAt(p).contentAt(c).texture(), std::string(TEMPO_PATH+"content_p") + std::to_string(p) + "_c" + std::to_string(c) + ".png");

//	std::string renderingDirectory = out_dir + "/" + name_noext + "_" + std::to_string(time((0)) + "/";
//	create_directory(renderingDirectory);
//	pProcessor.saveRenderingPack(renderingDirectory);

	ImageType output = pProcessor.generate().texture();
	IO::save01_in_u8(output, out_dir + "/output_ctexch_" + name_noext + "_" + std::to_string(time(nullptr)) + ".png");

	return 0;
}
