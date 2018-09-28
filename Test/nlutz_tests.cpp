#include <stdlib.h>
#include "ASTex/easy_io.h"
#include "ASTex/utils.h"
#include "ASTex/Stamping/stamper.h"
#include "ASTex/histogram.h"
#include "ASTex/rpn_utils.h"
#include "ASTex/mipmap.h"
#include "ASTex/texton_io.h"
#include "ASTex/imageviewer.h"
#include "ASTex/ContentExchange/patchProcessor.h"
#include "ASTex/ContentExchange/atlas.h"
#include "ASTex/PCTS/pcts.h"

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

    double sumW=0;
    double S=0;
    int max_dist=60.0, variance=2.0;
    if(argc>2)
        max_dist=atoi(argv[2]);
    if(argc>3)
        variance=atof(argv[3]);
    max_dist = max_dist%2 == 1 ? max_dist:max_dist-1; //must be odd
    ImageGrayd img, result;
    IO::loadu8_in_01(img, std::string(MY_PATH) + argv[1]);

    //building gaussian kernel of size max_dist
    ImageGrayd gaussianKernel;
    gaussianKernel.initItk(max_dist, max_dist, false);
    gaussianKernel.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
    {
        int middle=int(gaussianKernel.width())/2;
        if(x==middle && y==middle)
            pix=0;
        else
            pix = exp( -0.5*(pow(sqrt((x-middle)*(x-middle)+(y-middle)*(y-middle)) / variance, 2)) ) / (variance*sqrt(2*M_PI));
        sumW += pix;
    });
    //normalisation of the kernel
    gaussianKernel.for_all_pixels([&] (ImageGrayd::PixelType &pix)
    {
        pix /= sumW;
        std::cout << pix << std::endl;
    });

    HistogramGrayd h(img);

    double mean = h.mean();
    double sigma = h.variance();
    //mean and sigma of image =/= mean and sigma of zone

    //computing S
    img.for_all_pixels([&] (ImageGrayd::PixelType &pix)
    {
        S += pix*pix;
    });

    S /= img.width()*img.height();
    S -= mean*mean;
    S = std::sqrt(S);

    auto lmbd_gi = [&] (const ImageGrayd& texture, int x, int y) -> double
    {
        int dx, dy, x2, y2;
        double gi = 0, w;

        x += texture.width();
        y += texture.height();

//        double mean, sigma;
//        //computing local mean

//        for(dx=-gaussianKernel.width()/2; dx<=gaussianKernel.width()/2; ++dx)
//            for(dy=-gaussianKernel.height()/2; dy<=gaussianKernel.height()/2; ++dy)
//            {
//                x2 = (x + dx)%texture.width();
//                y2 = (y + dy)%texture.height();
//                mean += texture.pixelAbsolute(x2, y2);
//            }
//        mean /= gaussianKernel.width()*gaussianKernel.height();

//        //computing local sigma

//        for(dx=-gaussianKernel.width()/2; dx<=gaussianKernel.width()/2; ++dx)
//            for(dy=-gaussianKernel.height()/2; dy<=gaussianKernel.height()/2; ++dy)
//            {
//                x2 = (x + dx)%texture.width();
//                y2 = (y + dy)%texture.height();
//                sigma += (texture.pixelAbsolute(x2, y2)-mean) * (texture.pixelAbsolute(x2, y2)-mean);
//            }
//        sigma /= gaussianKernel.width()*gaussianKernel.height()-1;

        //computing Gi

        double denom_gi=0;

        for(dx=-gaussianKernel.width()/2; dx<=gaussianKernel.width()/2; ++dx)
            for(dy=-gaussianKernel.height()/2; dy<=gaussianKernel.height()/2; ++dy)
            {
                x2 = (x + dx)%texture.width();
                y2 = (y + dy)%texture.height();
                w = gaussianKernel.pixelAbsolute(dx + gaussianKernel.width()/2, dy + gaussianKernel.height()/2);
                gi += texture.pixelAbsolute(x2, y2) * w;

                denom_gi += w*w;
            }
        gi -= mean;
        denom_gi *= img.width()*img.height();
        denom_gi -= 1;
        denom_gi /= img.width()*img.height()-1;
        denom_gi = S*std::sqrt(denom_gi);

        gi /= denom_gi;
        return gi;

    };

    result.initItk(img.width(), img.height());
    result.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
    {
        pix = lmbd_gi(img, x, y);
        //std::cout << "gi(" << x << ", " << y << "): " << pix << std::endl;
        //pix+=1;
        //pix*=0.5;
        pix = std::abs(pix);
        pix = pix > 1.0 ? 1.0 : pix < 0 ? 0 : pix;
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

int test_contentFiltering(int argc, char **argv)
{
    if(argc < 3)
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " <in_texture> <in_patchMap>" << std::endl;
        return EXIT_FAILURE;
    }

    int k, xShift, yShift;
    ImageRGBd im_in, im_patches, im_out;

    std::vector<ImageRGBAd> patchesVector;
    std::vector<Eigen::Vector2f> translationsVector;

    ImageRGBAd im_fullMipmap;
    ImageRGBd im_groundTruthMipmap;
    ImageRGBAd im_globalMipmap;

    if(!IO::loadu8_in_01(im_in, std::string(MY_PATH)+argv[1]))
        return 1;
    if(!IO::loadu8_in_01(im_patches, std::string(MY_PATH)+argv[2]))
        return 1;

    HistogramRGBd patchesHisto(im_patches);
    patchesVector.reserve(patchesHisto.binsNumber());
    translationsVector.reserve(patchesHisto.binsNumber());

    Mipmap<ImageRGBd> groundTruthMipmap(im_in);
    groundTruthMipmap.setMode(ANISOTROPIC);
    groundTruthMipmap.generate();
    groundTruthMipmap.fullMipmap(im_groundTruthMipmap);
    im_globalMipmap.initItk(im_groundTruthMipmap.width(), im_groundTruthMipmap.height(), true);

    IO::save01_in_u8(im_groundTruthMipmap, std::string(MY_PATH)+"mipmap_full_groundTruth.png");

    k=0;

    for(HistogramRGBd::const_iterator c_it=patchesHisto.begin(); c_it!=patchesHisto.end(); ++c_it)
    {
        ImageRGBAd img;
        ImageRGBd::PixelType binValue=(*c_it).first;
        srand(k+42);
        xShift=rand();
        yShift=rand();
        patchesVector.push_back(img);
        patchesVector.back().initItk(im_in.width(), im_in.height(), true);
        patchesVector.back().for_all_pixels([&] (ImageRGBAd::PixelType &pix, int x, int y)
        {
            if( binValue == im_patches.pixelAbsolute(x, y) )
            {
                for(int i=0; i<3; ++i)
                    pix[i]=im_in.pixelAbsolute((x+xShift)%im_in.width(), (y+yShift)%im_in.height())[i];
                pix.SetAlpha(1.0);
            }
        });

        Mipmap<ImageRGBAd> mipmap(patchesVector.back());
        mipmap.setMode(ANISOTROPIC);
        mipmap.generate();
        mipmap.fullMipmap(im_fullMipmap);
        IO::save01_in_u8(im_fullMipmap, std::string(MY_PATH)+"mipmap_" + std::to_string(k++) + ".png");
        im_globalMipmap.for_all_pixels( [&] (ImageRGBAd::PixelType &pix, int x, int y)
        {
            pix += im_fullMipmap.pixelAbsolute(x, y);
        });
    }

    std::string filename_source=argv[1];
    std::string name_file = IO::remove_path(filename_source);
    std::string name_noext = IO::remove_ext(name_file);

    IO::save01_in_u8(im_globalMipmap, std::string(MY_PATH)+name_noext+"_mipmap_full_computed.png");

    return 0;
}

int test_contentExchangeFiltering(int argc, char **argv)
{
    if(argc < 3)
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " <in_texture> <in_patchMap>" << std::endl;
        return EXIT_FAILURE;
    }

    ImageRGBd im_in, im_patches;

    if(!IO::loadu8_in_01(im_in, std::string(MY_PATH)+argv[1]))
        return 1;
    if(!IO::loadu8_in_01(im_patches, std::string(MY_PATH)+argv[2]))
        return 1;

    std::string out_dir=std::string(MY_PATH)+"contentMipmaps/";
    create_directory(out_dir);

    ContentExchange::PatchProcessor<ImageRGBd> pProcessor(im_in);
    //pProcessor.setFilteringMode(ANISOTROPIC);
    pProcessor.setFilteringMode(NO_FILTER);
    //add some patches here
    pProcessor.initializePatchesFromImageRGB(im_patches);
    //patchProcessor.initializePatchesRegularGrid(64);
    pProcessor.initializeContents();
    //add some contents there
    pProcessor.debug_setRandomContents(3);
    std::cout << "memory cost of one pixel: " << std::to_string(sizeof(ImageRGBd::PixelType)) << std::endl;
    std::cout << "total memory cost: " << std::to_string(pProcessor.analysis_getGPUMemoryCost()) << std::endl;

    unsigned k=0;

    create_directory(out_dir + "/rendering");
    pProcessor.saveRenderingPack(out_dir + "/rendering");

    //TEST 1: EVERY CONTENTS OF EVERY PATCH, + EVERY ALPHA OF EVERY PATCH, + EVERY POSSIBLE COMBINATION PATCH+CONTENT

    if(pProcessor.filteringMode() == ANISOTROPIC)
    for(typename ContentExchange::PatchProcessor<ImageRGBd>::iterator it=pProcessor.begin(); it!=pProcessor.end(); ++it, ++k)
    {
        ImageRGBAd im_out;
        typename ContentExchange::Patch<ImageRGBd> &patch=(*it);
        for(unsigned i=0; i<pProcessor.numberMipmapsWidth(); ++i)
            for(unsigned j=0; j<pProcessor.numberMipmapsHeight(); ++j)
            {
                IO::save01_in_u8(patch.mipmap(i, j), out_dir + "alpha_p" + std::to_string(k) + "_c0"
                                 + "_mw" + std::to_string(i) + "_mh" + std::to_string(j) + ".png");
                for(unsigned l=0; l<patch.nbContents(); ++l)
                {
                    const ImageRGBd &contentMipmap=patch.contentAt(l).mipmap(i, j);
                    IO::save01_in_u8(contentMipmap, out_dir + "content_p" + std::to_string(k) + "_c" + std::to_string(l)
                                     + "_mw" + std::to_string(i) + "_mh" + std::to_string(j) + ".png");
                    im_out.initItk(contentMipmap.width(), contentMipmap.height(), true);
                    contentMipmap.for_all_pixels([&] (const ImageRGBd::PixelType &pix, int x, int y)
                    {
                        ImageRGBAd::PixelType p;
                        p.SetRed(pix.GetRed());
                        p.SetGreen(pix.GetGreen());
                        p.SetBlue(pix.GetBlue());
                        p.SetAlpha(patch.mipmap(i, j).pixelAbsolute(x, y));
                        im_out.pixelAbsolute(x, y)=p;
                    });
                    IO::save01_in_u8(im_out, out_dir + "mipmap_p" + std::to_string(k) + "_c" + std::to_string(l)
                                     + "_mw" + std::to_string(i) + "_mh" + std::to_string(j) + ".png");
                }
            }
    }

    //TEST 2 : OUTPUT MIPMAP WITH EVERY DEFAULT CONTENTS

    ImageRGBd im_out;
    Mipmap<ImageRGBd> outputMipmap(pProcessor.generate(512, 512));

    size_t s=0;
    for(unsigned m=0; m<outputMipmap.numberMipmapsWidth(); ++m)
    {
        if(outputMipmap.mode()==ISOTROPIC)
            s+=outputMipmap.mipmap(m, m).width()*outputMipmap.mipmap(m, m).height()*sizeof(ImageRGBd::PixelType);
        else
            for(unsigned l=0; l<outputMipmap.numberMipmapsHeight(); ++l)
            {
                int access;
                access = pProcessor.analysis_getNumberOfTextureAccessForMipmap(m, l);
                std::cout << "Mipmap(" << m << ", " << l << "): ";
                std::cout << "The number of texture access for computing this mipmap is " << std::to_string(access)
                          << ", which means there are " << std::to_string(access/(double(outputMipmap.mipmap(m, l).width()*outputMipmap.mipmap(m, l).height())))
                          << " reads per pixel, in mean." << std::endl;
                s+=outputMipmap.mipmap(m, l).width()*outputMipmap.mipmap(m, l).height()*sizeof(ImageRGBd::PixelType);
            }
    }

    std::cout << "Classic mipmap's memory cost: " << s << std::endl;

    outputMipmap.fullMipmap(im_out);
    IO::save01_in_u8(im_out, std::string(MY_PATH)+"contentMipmaps/_fullMipmap.png");
    outputMipmap.setTexture(im_in);
    outputMipmap.generate();
    outputMipmap.fullMipmap(im_out);
    IO::save01_in_u8(im_out, std::string(MY_PATH)+"contentMipmaps/_fullMipmapOfInput.png");

    ContentExchange::Atlas<ImageRGBd> atlas(pProcessor);
    atlas.generate(1);
    im_out = atlas.generatedImage();

    IO::save01_in_u8(im_out, std::string(MY_PATH) + "contentMipmaps/_atlasOffline.png");

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

char *name = "bricks";
char fname[256];

int test_pcts(int argc, char **argv)
{
    std::cout << "PCTS started" << std::endl;

    ImageRGBd image, guid, seg;
	ImageRGBd mask, Ipos, Ineg, I2pos, I2neg;
	ImageRGBd synth, binary, stencil;
	//set this image and macro PCTS_DEBUG_DIRECTORY in pcts.h
	sprintf(fname, "E:/developpement/AsTex/AsTex/Data/%s.png", name);
    IO::loadu8_in_01(image, fname);
	sprintf(fname, "E:/developpement/AsTex/AsTex/Data/%s_mask.png", name);
	IO::loadu8_in_01(mask, fname);
	sprintf(fname, "E:/developpement/AsTex/AsTex/Data/%s_init_Binary_warped_specific_DT.png", name);
	IO::loadu8_in_01(Ipos, fname);
	sprintf(fname, "E:/developpement/AsTex/AsTex/Data/%s_init_Binary_warped_specific_DT_neg.png", name);
	IO::loadu8_in_01(Ineg, fname);
	sprintf(fname, "E:/developpement/AsTex/AsTex/Data/%s_rigidity_RGB.png", name);
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
	sprintf(fname, "E:/developpement/AsTex/AsTex/Data/%s_Binary_warped_specific_DT.png", name);
	IO::loadu8_in_01(I2pos, fname);
	sprintf(fname, "E:/developpement/AsTex/AsTex/Data/%s_Binary_warped_specific_DT_neg.png", name);
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
	sprintf(fname, "E:/developpement/AsTex/AsTex/Data/%s_C_10.png", name);
	IO::loadu8_in_01(synth, fname);
	sprintf(fname, "E:/developpement/AsTex/AsTex/Data/%s_C_10_bin.png", name);
	IO::loadu8_in_01(binary, fname);

	ASTex::Pcts<ImageRGBd> pcts;
    pcts.setTexture(image);
    //pcts.setWidth(800);
    //pcts.setHeight(800);
    pcts.setNbSamplesNNM(20);
    pcts.setNbRefinementsNNM(2);
    pcts.setRadiusScaleNNM(15);
	pcts.setLabel(mask, 0.5);
	pcts.setGuidance(guid, seg, 0.95, 0.01);
	pcts.setMask(synth, binary);
	pcts.setStencil(stencil, 1.0); // weight between 0 and 1
    IO::save01_in_u8(pcts.generate(), "E:/developpement/AsTex/AsTex/Data/bricks_pcts.png");

    std::cout << "PCTS ended" << std::endl;
    return 0;
}

int main( int argc, char **argv )
{
    //return test_genet(argc, argv);
    //return test_texton(argc, argv);
    //return test_autoconvolutionSpectrum(argc, argv);
    //return test_contentExchangeFiltering(argc, argv);
    //return test_wendling(argc, argv);
    //return test_getis_gi(argc, argv);
    //return test_easy(argc, argv);
    return test_pcts(argc, argv);
}
