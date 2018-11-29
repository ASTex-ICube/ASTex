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
        ImageRGBu8 im_in;
        im_in.load(std::string(std::string(MY_PATH)+argv[i]));

        std::string out_dir=std::string(MY_PATH)+"contentExchangeRendering";
        std::string name_file = IO::remove_path(std::string(MY_PATH)+argv[i]);
        std::string name_noext = IO::remove_ext(name_file);
        create_directory(out_dir);

        ContentExchange::PatchProcessor<ImageRGBu8> pProcessor(im_in);
        pProcessor.setFilteringMode(ISOTROPIC);
        pProcessor.setNbContentsPerPatch(10);
        pProcessor.fullProcess_oldMethod();

        std::string renderingDirectory = out_dir + "/" + std::to_string(std::time(0)) + "_" + name_noext + "/";
        create_directory(renderingDirectory);
        pProcessor.saveRenderingPack(renderingDirectory);
    }

    return 0;
}

/**
 * @brief test_texton is a template for making texton noise.
 * Change the variable MY_PATH to the patch you have your images.
 * You can safely change the following lines:
 * sampler.setNbPoints(300);
 *      ^ which specifies the number of times texton noise hits the output texture
 * stamp.setInterpolationRule(Stamping::StampDiscrete<ImageRGBd>::BILINEAR);
 *      ^ which specifies the interpolation rule (we shoot between pixels)
 * tamponneur.setPeriodicity(false);
 *      ^ which specifies whether the output image is allowed to be periodic or not
 * int W=1024, H=1024;
 *      ^ which specifies the output image width and height.
 * @param argc main argc
 * @param argv main argv
 * @return 0 (always !)
 */
int test_texton(int argc, char **argv)
{
    if( argc < 3 )
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << "<input image> <texton image/file>" << std::endl;

        return EXIT_FAILURE;
    }

#ifdef USE_QWIDGETS
    QApplication app(argc, argv);
    std::setlocale(LC_ALL,"C");
#endif

    std::string filename_source=argv[1];
    std::string name_file = IO::remove_path(filename_source);
    std::string name_noext = IO::remove_ext(name_file);

    std::string input_noext=name_noext;

    filename_source=argv[2];
    name_file = IO::remove_path(filename_source);
    name_noext = IO::remove_ext(name_file);

    //Loading sample
    //Loading im_in

    ImageRGBd im_in, im_out, im_rpn, im_sample, im_texton;
    IO::loadu8_in_01(im_in, std::string(MY_PATH)+argv[1]);

    //warning; these import functions turn a texton file into a vizualizable image
    if(!import_texton(im_texton, MY_PATH+filename_source))
    {
        assert(import_texton_from_png(im_texton, MY_PATH+filename_source));
    }

    HistogramRGBd histo_texton(im_texton);
    ImageRGBd::PixelType mean;
    for(int i=0; i<3; ++i)
        mean[i] = histo_texton.mean(i);
    //transformation image -> texton
    im_texton.for_all_pixels([&] (ImageRGBd::PixelType &pix) {
        pix -= mean;
        pix = pix * (1.0/std::sqrt(im_texton.width()*im_texton.height()));
    });

    //Testing

    Stamping::SamplerUniform sampler;
    sampler.setNbPoints(400); //< you can change that
    Stamping::StampDiscrete<ImageRGBd> stamp(im_texton);
    stamp.setInterpolationRule(Stamping::StampDiscrete<ImageRGBd>::BILINEAR); //< you can change that too

    Stamping::StamperTexton<ImageRGBd> tamponneur(&sampler, &stamp);

    tamponneur.setPeriodicity(false);
    tamponneur.setUseMargins(true);

    int W=512, H=512;

    im_out = tamponneur.generate(W, H);

    //transformation texton -> image
    im_out.for_all_pixels([&] (ImageRGBd::PixelType &pix) {
        //pix = pix * std::sqrt(im_texton.width()*im_texton.height()); //mistake: no need to re-normalize on top
        pix += mean;
    });

    im_sample.initItk(W, H, true);


    std::vector<Eigen::Vector2f> pointArray = sampler.generate();
    for(std::vector<Eigen::Vector2f>::iterator it=pointArray.begin(); it!=pointArray.end(); ++it)
    {
        int i = im_sample.width() * (*it)[0]; //i & j: single point coordinates in im_out
        int j = im_sample.height() * (*it)[1];

        im_sample.pixelAbsolute(i, j)=1.0;
    }

    //colored_RPN(im_in, im_rpn, RGB_SPACE, NORMAL, 0, 0, false, true, false, 1.0);

    auto clamp = [&] (ImageRGBd::PixelType &pix) { for(int i=0; i<3; ++i) pix[i] = pix[i] > 1.0 ? 1.0 : (pix[i] < 0.0 ? 0.0 : pix[i]); };
    im_out.for_all_pixels(clamp);
    //im_rpn.for_all_pixels(clamp);

    HistogramRGBd histo_in(im_in);
    HistogramRGBd histo_out(im_out);
    //HistogramRGBd histo_rpn(im_rpn);

    std::cout << "Mean of input is: (" << std::to_string(histo_in.mean(0)) << ", " << std::to_string(histo_in.mean(1)) << ", " << std::to_string(histo_in.mean(2)) << ")" << std::endl;
    std::cout << "Variance of input is: (" << std::to_string(histo_in.covariance(0, 0)) << ", " << std::to_string(histo_in.covariance(1, 1)) << ", " << std::to_string(histo_in.covariance(2, 2)) << ")" << std::endl;

    std::cout << "==================" << std::endl;

    std::cout << "Mean of output is: (" << std::to_string(histo_out.mean(0)) << ", " << std::to_string(histo_out.mean(1)) << ", " << std::to_string(histo_out.mean(2)) << ")" << std::endl;
    std::cout << "Variance of output is: (" << std::to_string(histo_out.covariance(0, 0)) << ", " << std::to_string(histo_out.covariance(1, 1)) << ", " << std::to_string(histo_out.covariance(2, 2)) << ")" << std::endl;

    std::cout << "==================" << std::endl;

//    std::cout << "Mean of rpn is: (" << std::to_string(histo_rpn.mean(0)) << ", " << std::to_string(histo_rpn.mean(1)) << ", " << std::to_string(histo_rpn.mean(2)) << ")" << std::endl;
//    std::cout << "Variance of rpn is: (" << std::to_string(histo_rpn.covariance(0, 0)) << ", " << std::to_string(histo_rpn.covariance(1, 1)) << ", " << std::to_string(histo_rpn.covariance(2, 2)) << ")" << std::endl;

//    std::cout << "==================" << std::endl;

//    std::cout << "Mean of texton is: (" << std::to_string(histo_texton.mean(0)) << ", " << std::to_string(histo_texton.mean(1)) << ", " << std::to_string(histo_texton.mean(2)) << ")" << std::endl;
//    std::cout << "Variance of texton is: (" << std::to_string(histo_texton.covariance(0, 0)) << ", " << std::to_string(histo_texton.covariance(1, 1)) << ", " << std::to_string(histo_texton.covariance(2, 2)) << ")" << std::endl;


    double zero[3];
    double one[3];
    for(int i=0; i<3; ++i) {zero[i]=0.0; one[i]=1.0;}

//    HistogramRGBBase<int> quantizedHisto;

//    quantizedHisto=histo_out.quantize(zero, one, 24);
//    quantizedHisto.saveHistogram(std::string(MY_PATH) + "out_" + input_noext + "_tn" + ".csv", 24);

//    quantizedHisto=histo_rpn.quantize(zero, one, 24);
//    quantizedHisto.saveHistogram(std::string(MY_PATH) + "out_" + input_noext + "_rpn" + ".csv", 24);

//    if(argc > 3) //add an original texton noise image to compare the results with it
//    {
//        filename_source=argv[3];
//        ImageRGBd im_tn;
//        IO::loadu8_in_01(im_tn, MY_PATH + filename_source);

//        HistogramRGBd histo_tn(im_tn);

//        std::cout << "==================" << std::endl;

//        std::cout << "Mean of provided texton noise is: (" << std::to_string(histo_tn.mean(0)) << ", " << std::to_string(histo_tn.mean(1)) << ", " << std::to_string(histo_tn.mean(2)) << ")" << std::endl;
//        std::cout << "Variance of provided texton noise is: (" << std::to_string(histo_tn.covariance(0, 0)) << ", " << std::to_string(histo_tn.covariance(1, 1)) << ", " << std::to_string(histo_tn.covariance(2, 2)) << ")" << std::endl;

//        quantizedHisto=histo_tn.quantize(zero, one, 24);
//        quantizedHisto.saveHistogram(std::string(MY_PATH) + "out_" + input_noext + "_tn_galerne" + ".csv", 24);
//    }

    //Saving sample
    //Saving im_out
    IO::save01_in_u8(im_sample, std::string(MY_PATH) + name_noext + "_sample.png");
    IO::save01_in_u8(im_out, std::string(MY_PATH) + name_noext + "_out_tn.png");
//    IO::save01_in_u8(im_rpn, std::string(MY_PATH) + input_noext + "_out_rpn.png");
    IO::save01_in_u8(im_texton, std::string(MY_PATH) + name_noext + "_texton.png");

#ifdef USE_QWIDGETS
    ImageViewer imgv_in("Source", &app, 0);
    imgv_in.setWindowTitle("Source");
    imgv_in.set_rgb01(im_in.getDataPtr(), im_in.width(),im_in.height(),1);
    imgv_in.show();

    ImageViewer imgv_sample("Sampling", &app, 1);
    imgv_sample.setWindowTitle("Sampling");
    imgv_sample.set_rgb01(im_sample.getDataPtr(), im_sample.width(), im_sample.height(),1);
    imgv_sample.show();

//    ImageViewer imgv_rpn("RPN", &app, 2);
//    imgv_rpn.setWindowTitle("RPN");
//    imgv_rpn.set_rgb01(im_rpn.getDataPtr(), im_rpn.width(),im_rpn.height(),1);
//    imgv_rpn.show();

    ImageViewer imgv_out("Texton noise", &app, 3);
    imgv_out.setWindowTitle("Texton noise");
    imgv_out.set_rgb01(im_out.getDataPtr(), im_out.width(), im_out.height(),1);
    imgv_out.show();

    return app.exec();
#else
    return 0;
#endif
}

int test_easy(int argc, char **argv)
{
    std::cout << sizeof(ImageRGBu8::DataType) << std::endl;
    ImageRGBu8::PixelType pix;
    pix[0]=17;
    pix[1]=17;
    pix[2]=17;
    std::cout << int(*((unsigned char *)&pix+1)) << std::endl;
    return 0;
}

char fname[256];

int test_pcts(int argc, char **argv)
{
    std::cout << "PCTS started" << std::endl;

    ImageRGBd image, guid, seg;
	ImageRGBd mask, Ipos, Ineg, I2pos, I2neg;
	ImageRGBd synth, binary, stencil;

    std::string input_directory = std::string(argv[1]);
    char *name = argv[2];

	//set this image and macro PCTS_DEBUG_DIRECTORY in pcts.h
    sprintf(fname, "%s/%s.png", input_directory.c_str(), name);
    IO::loadu8_in_01(image, fname);
    sprintf(fname, "%s/mask.png", input_directory.c_str());
	IO::loadu8_in_01(mask, fname);
    sprintf(fname, "%s/init_Binary_warped_specific_DT.png", input_directory.c_str());
    IO::loadu8_in_01(Ipos, fname);
    sprintf(fname, "%s/init_Binary_warped_specific_DT_neg.png", input_directory.c_str());
	IO::loadu8_in_01(Ineg, fname);
    sprintf(fname, "%s/rigidity_RGB.png", input_directory.c_str());
    IO::loadu8_in_01(stencil, fname);
	seg.initItk(Ipos.width(), Ipos.height());
	seg.for_all_pixels([&](ImageRGBd::PixelType &pix, int x, int y)
	{
		double col[3];
		col[0] = pow(Ipos.pixelAbsolute(x, y)[0],0.25);
		col[1] = pow(Ineg.pixelAbsolute(x, y)[0],0.25);
		col[2] = 0.0;
		pix = ImageRGBd::PixelType(col);
	});
    sprintf(fname, "%s/Binary_warped_specific_DT.png", input_directory.c_str());
	IO::loadu8_in_01(I2pos, fname);
    sprintf(fname, "%s/Binary_warped_specific_DT_neg.png", input_directory.c_str());
	IO::loadu8_in_01(I2neg, fname);
	guid.initItk(I2pos.width(), I2pos.height());
	guid.for_all_pixels([&](ImageRGBd::PixelType &pix, int x, int y)
	{
		double col[3];
		col[0] = pow(I2pos.pixelAbsolute(x, y)[0],0.25);
		col[1] = pow(I2neg.pixelAbsolute(x, y)[0],0.25);
		col[2] = 0.0;
		pix = ImageRGBd::PixelType(col);
	});
    sprintf(fname, "%s/%s_C_10.png", input_directory.c_str(), name);
	IO::loadu8_in_01(synth, fname);
    sprintf(fname, "%s/%s_C_10_bin.png", input_directory.c_str(), name);
	IO::loadu8_in_01(binary, fname);

	ASTex::Pcts<ImageRGBd> pcts;
    pcts.setTexture(image);
    //pcts.setWidth(800);
    //pcts.setHeight(800);


    unsigned    seed=0,
                bs = 4,
                sl = 5,
                samples = 60,
                ref = 1,
                radius = 30;

     double     labelW = 0.9,
                guidanceW = 0.8,
                stencilW = 1.0;

    srand(seed);
    pcts.setLvl0BlockSize(bs);
    pcts.setMinimumSizeLog(sl);
    pcts.setNbSamplesNNM(samples);
    pcts.setNbRefinementsNNM(ref);
    pcts.setRadiusScaleNNM(radius);
    pcts.setLabel(mask, labelW);
    pcts.setGuidance(guid, seg, guidanceW, 1.0);
	pcts.setMask(synth, binary);
    pcts.setStencil(stencil, stencilW); // weight between 0 and 1
    IO::save01_in_u8(pcts.generate(), input_directory + "/out_pcts_seed" + std::to_string(seed)
                     + "_bs" + std::to_string(bs)
                     + "_sl" + std::to_string(sl)
                     + "_samples" + std::to_string(samples)
                     + "_ref" + std::to_string(ref)
                     + "_labelW" + std::to_string(labelW)
                     + "_guidanceW" + std::to_string(guidanceW)
                     + "_stencilW" + std::to_string(stencilW)
                     + ".png");

    std::cout << "PCTS ended" << std::endl;
    return 0;
}

int test_benchmarkingContentExchange(int argc, char **argv)
{
    for(int i=1; i<argc; ++i)
    {
        ImageRGBu8 input;
        input.load(std::string("/home/nlutz/img/") + argv[i]);
        std::string name_file = IO::remove_path(argv[i]);
        std::string name_noext = IO::remove_ext(name_file);
        std::string pcts_file = MY_PATH + "pcts_" + name_noext + ".txt";
        ContentExchangeBenchmarker ceb;
        ceb.setInput(input);

        ceb.setNbOutputs(3);
        ceb.setOutputSize(1024, 1024);
        ceb.setRoot(std::string("/home/nlutz/ieee2019/") + name_noext);
        ceb.setOutputDirectories("usingOldMethod", "usingRandomPatchesAndContents", "usingGetisGI", "usingPCTS");
        ceb.setGenerateOld(true);
        ceb.setGenerateRandom(false);
        ceb.setGeneratePCTS(false);
        ceb.setGenerateGetisGI(false);
        ceb.setPCTSArgumentsFilePath(pcts_file);

        ceb.generate();
    }
    return 0;
}

int main( int argc, char **argv )
{
    //return test_genet(argc, argv);
    //return test_texton(argc, argv);
    //return test_autoconvolutionSpectrum(argc, argv);
    return test_contentExchangeRenderingPack(argc, argv);
    //return test_wendling(argc, argv);
    //return test_getis_gi(argc, argv);
    //return test_easy(argc, argv);
    //return test_pcts(argc, argv);
    //return test_benchmarkingContentExchange(argc, argv);
}
