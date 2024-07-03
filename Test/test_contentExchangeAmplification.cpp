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

	//nb de patch
	//taille de la tuile
	//largeur du blend
	//quels exemples utilisés + transformations
	//résolution de la guidance

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

	auto lmbd_centerGradient = [] (ImageType &gradientImage)
	{
		gradientImage.for_all_pixels([&] (ImageType::PixelType &pix)
		{
			pix[0] -= 0.5;
			pix[1] -= 0.5;
		});
	};
	auto lmbd_uncenterGradient = [] (ImageType &gradientImage)
	{
		gradientImage.for_all_pixels([&] (ImageType::PixelType &pix)
		{
			pix[0] += 0.5;
			pix[1] += 0.5;
		});
	};

	auto lmbd_normalizeImage = [] (ImageType &gradientImage)
	{	//I use this to normalize between 0 and 1 even when the image is 16 bit
		//but feel free to come up with something easier and more robust
		bool is16bit = false;
		gradientImage.for_all_pixels([&] (ImageType::PixelType &pix)
		{
			if(pix[0]>255 || pix[1]>255)
				is16bit=true;
		});
		gradientImage.for_all_pixels([&] (ImageType::PixelType &pix)
		{
			for(unsigned i=0; i<3; ++i)
			{
				if(is16bit)
					pix[i]/=65535.0;
				else
					pix[i]/=255.0;
			}
		});
	};

	ImageType im_tile;
	im_tile.load(argv[2]);
	lmbd_normalizeImage(im_tile);

	ImageRGBd im_guidance;
	im_guidance.load(argv[3]);
	lmbd_normalizeImage(im_guidance);

	std::string out_dir = argv[1];
	std::string name_file = IO::remove_path(argv[3]);
	std::string name_noext = IO::remove_ext(name_file);
	create_directory(out_dir);

	///Argument file
	std::ifstream ifs_args(argv[4]);
	if(!ifs_args)
	{
		std::cerr << "Error: could not open argument file" << std::endl;
		exit(EXIT_FAILURE);
	}
	unsigned nb_patches, nb_contents, output_width, output_height, seed, gradientManipulation;
	ifs_args >> nb_patches >> nb_contents >> output_width >> output_height >> seed >> gradientManipulation;

	ASTex::TexturePool<ImageType, double> texturePool;

	ContentExchange::MultiTextureProcessor<ImageType> pProcessor(im_tile);
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
	//texturePool.addTexture(im_tile, true);
	texturePool.setTexturesAreGradient(gradientManipulation);
//	texturePool.addFunction([](ImageType &image)
//	{//this one makes the gradient scale two times lower
//		image.for_all_pixels([&] (ImageType::PixelType &pix)
//		{
//			for(unsigned i=0; i<2; ++i)
//			{
//				if(pix[i]<0.5)
//				{
//					pix[i]=0.5-(0.5-pix[i])*0.5;
//				}
//				else
//				{
//					pix[i]=0.5+(pix[i]-0.5)*0.5;
//				}
//			}
//		});
//	});

	if(gradientManipulation)
	{
		lmbd_centerGradient(im_tile);
		lmbd_centerGradient(im_guidance);
	}

	for(int i=5; i<argc - (argc%2==0 ? 2 : 0); i+=2)
	{
		ImageType im_texture;

		im_texture.load(argv[i]);
		lmbd_normalizeImage(im_texture);
		int periodicity = atoi(argv[i+1]);
		lmbd_centerGradient(im_texture);
		texturePool.addTexture(im_texture, bool(periodicity));
	}
	if(argc%2==0)
	{
		ImageType im_texture;
		im_texture.load(argv[argc-1]);
		lmbd_normalizeImage(im_texture);
		lmbd_centerGradient(im_texture);
		texturePool.addTexture(im_texture, false);
	}
	texturePool.generate();
	for(unsigned i=0; i<texturePool.size(); ++i)
	{
		ImageRGBd texture;
		texture.initItk(texturePool[i].texture.width(), texturePool[i].texture.height());
		texture.copy_pixels(texturePool[i].texture);
		lmbd_uncenterGradient(texture);
		IO::save01_in_u8(texturePool[i].texture, TEMPO_PATH+"expTexture.png");
		IO::save01_in_u8(texture, TEMPO_PATH+"expTexture_i" + std::to_string(i) + ".png");
//		IO::save01_in_u8(texturePool[i].boundaries, TEMPO_PATH+"inBoundsMipmap_i" + std::to_string(i) + ".png");
	}

	pProcessor.setSeed(seed == 0 ? time(nullptr) : seed);

	pProcessor.setFilteringMode(ISOTROPIC);
	pProcessor.setNbContentsPerPatch(nb_contents);
	pProcessor.setOutputSize(output_width == 0 ? unsigned(im_guidance.width()) : output_width,
							 output_height == 0 ? unsigned(im_guidance.height()) : output_height);
	pProcessor.patches_initRandom(nb_patches);
	pProcessor.patches_dilate(false);
	pProcessor.patches_dilate(false);
	pProcessor.contents_initDefault();
	pProcessor.setTexturePool(texturePool);

	Mipmap<ImageType> mipmapOutput = pProcessor.generateFromExemplar(im_guidance);
	ImageRGBd outputTexture = mipmapOutput.texture();
	lmbd_uncenterGradient(outputTexture);
	outputTexture.for_all_pixels([&] (ImageType::PixelType &pix)
	{
		for(unsigned i=0; i<3; ++i)
			pix[i] = std::min(std::max(0.0, pix[i]), 1.0);
	});
//	IO::save01_in_u8(outputTexture, out_dir + "/output_ctexchAmp_projection_" + name_noext + std::to_string(time(nullptr)) + ".png");
	ImageRGBu16 output16;
	output16.initItk(outputTexture.width(), outputTexture.height());
	output16.for_all_pixels([&] (ImageRGBu16::PixelType &pix, int x, int y)
	{
		for(unsigned i=0; i<2; ++i)
			pix[i] = uint16_t(outputTexture.pixelAbsolute(x, y)[i]*65535);
		pix[2] = 65535;
	});
	output16.save(out_dir + "/output_ctexchAmp_projection_" + name_noext + "_tw" + std::to_string(im_tile.width()) + "_w" + std::to_string(output16.width()) + "_" + std::to_string(time(nullptr)) + ".png");

	return 0;
}
