#include <stdlib.h>
#include "ASTex/easy_io.h"
#include "ASTex/utils.h"
#include "ASTex/Stamping/stamper.h"
#include "ASTex/histogram.h"
#include "ASTex/rpn_utils.h"
#include "ASTex/mipmap.h"
#include "ASTex/texton_io.h"
#include "ASTex/imageviewer.h"
#include "ASTex/ContentExchange/atlas.h"
#include "ASTex/PCTS/pcts.h"
#include "ASTex/filters.h"
#include "ASTex/ContentExchange/benchmarker.h"

// v define your own
#define MY_PATH std::string("/home/nlutz/img/")

using namespace ASTex;

int test_getis_gi(int argc, char **argv)
{
	if(argc < 2)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " <in_texture> [max_dist (pixels)] [variance]" << std::endl;
		return EXIT_FAILURE;
	}

	ImageGrayd img;
	IO::loadu8_in_01(img, std::string(MY_PATH) + argv[1]);
	int max_dist = argc > 2 ? atoi(argv[2]) : 60;
	double variance = argc > 3 ? atof(argv[3]) : 1.0;

	FilterGetisGI<ImageGrayd>::Pointer fggi = FilterGetisGI<ImageGrayd>::New();
	fggi->setKernelSize(max_dist);
	fggi->setKernelVariance(variance);
	fggi->SetInput(img.itk());
	fggi->Update();
	ImageGrayd result(fggi->GetOutput());

	result.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
	{
		pix = std::max(0.0, std::min(1.0, (std::abs(pix) * 0.1)));
	});

	std::string filename_source=argv[1];
	std::string name_file = IO::remove_path(filename_source);
	std::string name_noext = IO::remove_ext(name_file);
	IO::save01_in_u8(result, std::string(MY_PATH)+ "gi_translated_" + name_noext + "_" + std::to_string(max_dist) + "_" + std::to_string(variance) + ".png");

	return 0;
}

int test_wendling(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " <in_texture> <in_patchMap>" << std::endl;
		return EXIT_FAILURE;
	}
	auto getVariance = [&] (const ImageGrayd& img){
			double variance =0;
			double mean = 0;

			img.for_all_pixels([&](const ImageGrayd::PixelType& p){
			   mean+=p;
			});

			mean /= img.width()*img.height();

			img.for_all_pixels([&](const ImageGrayd::PixelType& p){
			   double diff_mean = mean-p;
			   variance += diff_mean*diff_mean;
			});
			return variance/(img.width()*img.height());
		};

		ImageGrayd img, result;
		IO::loadu8_in_01(img, std::string(MY_PATH) + argv[1]);

		double mean = 0;
		ImageGrayd tmp(img.width(),img.height());

		//initialisation de l'image de travail (copie de l'image de base
		tmp.for_all_pixels([&](ImageGrayd::PixelType& p,int x,int y){
			p = img.pixelAbsolute(x,y);

			//Calcul de la moyenne de l'image
			mean+=p;
		});

		mean /= img.width()*img.height();

		//val-moyenne pour le calcul de l'autocorrelation
		tmp.for_all_pixels([&](ImageGrayd::PixelType& p){
			p-=mean;
		});

		//Calcul de l'autocorrelation
		result.initItk(tmp.width(), tmp.height());
		Fourier::autoCorrelation_full_size(tmp,result);

		double var_img = getVariance(img);

		//Normalisation par la variance de l'image pour val dans [-1;1]
		result.for_all_pixels([&](ImageGrayd::PixelType& p){
			p/=var_img;
		});

		return 0;
}

int test_contentExchangeRenderingPack(int argc, char **argv)
{
    if(argc < 2)
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " <in_texture>" << std::endl;
        return EXIT_FAILURE;
    }

    for(int i=1; i<argc; ++i)
    {
        srand(0);
		ImageRGBf im_in;
		//im_in.load(std::string(std::string(MY_PATH)+argv[i]));
		IO::loadu8_in_01(im_in, std::string(std::string(MY_PATH)+argv[i]));

        std::string out_dir=std::string(MY_PATH)+"contentExchangeRendering";
        std::string name_file = IO::remove_path(std::string(MY_PATH)+argv[i]);
        std::string name_noext = IO::remove_ext(name_file);
        create_directory(out_dir);

		ContentExchange::PatchProcessor<ImageRGBf> pProcessor(im_in);
        pProcessor.setFilteringMode(ISOTROPIC);
        pProcessor.setNbContentsPerPatch(6);
        //pProcessor.fullProcess_oldMethod(); //for nerds
        pProcessor.patches_initRandom(32);
        pProcessor.patches_dilate(false);
        pProcessor.contents_initDefault();
        pProcessor.contents_initRandom();

//        pProcessor.patchAt(1).contentAt(0).mipmap(0, 0).save("/home/nlutz/mipmap1.png");
//        pProcessor.patchAt(1).contentAt(0).mipmap(1, 1).save("/home/nlutz/mipmap2.png");
//        pProcessor.patchAt(1).contentAt(0).mipmap(2, 2).save("/home/nlutz/mipmap3.png");
//        pProcessor.patchAt(1).contentAt(0).mipmap(3, 3).save("/home/nlutz/mipmap4.png");

		std::string renderingDirectory = out_dir + "/" + std::to_string(std::time(0)) + "_random_" + name_noext + "/";

		ImageRGBf im_patches;
        im_patches.initItk(im_in.width(), im_in.height());
		im_patches.for_all_pixels([&] (ImageRGBf::PixelType &pix, int x, int y)
        {
            ContentExchange::word64 w = pProcessor.patchMapMipmap().texture().pixelAbsolute(x, y);
            pix = std::ceil(log2(w)) != std::floor(log2(w)) ? 255 : 0;
        });

		im_patches.save("/home/nlutz/im_seams.png");

		create_directory(renderingDirectory);
		pProcessor.saveRenderingPack(renderingDirectory, false);
	}

	return 0;
}

int test_easy(int argc, char **argv)
{
	ImageRGBu8 input;
	input.load(argv[1]);

	std::cout << bilinear_interpolation(input, 0.5, 0.5, true) << std::endl;

	return 0;
}

int test_benchmarkingContentExchange(int argc, char **argv)
{
    for(int i=1; i<argc; ++i)
    {
        ImageRGBu8 input;
		input.load(std::string(argv[i]));
        std::string name_file = IO::remove_path(argv[i]);
        std::string name_noext = IO::remove_ext(name_file);
        std::string pcts_file = MY_PATH + "pcts_" + name_noext + ".txt";
        ContentExchangeBenchmarker ceb;
		srand(42);
        ceb.setInput(input);
		ceb.setNbContentsPerPatch(6);
        ceb.setNbOutputs(3);
		ceb.setOutputSize(512, 512);
        ceb.setRoot(std::string("/home/nlutz/eg2019_contentExchange/") + name_noext);
        ceb.setOutputDirectories("usingOldMethod", "usingRandomPatchesAndContents", "usingGetisGI", "usingPCTS");
		ceb.setGenerateOld(true);
		ceb.setGenerateRandom(true);
        ceb.setGeneratePCTS(false);
        ceb.setGenerateGetisGI(false);
		//ceb.setPCTSArgumentsFilePath(pcts_file);

        ceb.generate();
    }
    return 0;
}

int main( int argc, char **argv )
{
    //return test_genet(argc, argv);
    //return test_texton(argc, argv);
    //return test_autoconvolutionSpectrum(argc, argv);
	//return test_contentExchangeRenderingPack(argc, argv);
    //return test_wendling(argc, argv);
    //return test_getis_gi(argc, argv);
	return test_easy(argc, argv);
	//return test_pcts(argc, argv);
	//return test_benchmarkingContentExchange(argc, argv);
}
