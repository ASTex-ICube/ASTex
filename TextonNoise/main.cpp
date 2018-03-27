#include <stdlib.h>
#include "ASTex/easy_io.h"
#include "ASTex/utils.h"
#include "Stamping/stamper.h"
#include "histogram.h"
#include "rpn_utils.h"
#include "mipmap.h"
#include "texton_io.h"

#include "imageviewer.h"

// v define your own
#define MY_PATH std::string("/home/nlutz/img/")

using namespace ASTex;

//int test_genet(int argc, char **argv)
//{
//    QApplication app(argc, argv);
//    std::setlocale(LC_ALL,"C");

//    std::string filename_source=argv[1];
//    std::string name_file = IO::remove_path(filename_source);
//    std::string name_noext = IO::remove_ext(name_file);

//    //std::string out_path = argv[1];

//    std::string input_noext=name_noext;

//    filename_source=argv[2];
//    name_file = IO::remove_path(filename_source);
//    name_noext = IO::remove_ext(name_file);

//    //Loading sample
//    //Loading im_in

//    ImageGrayd im_in, im_out;
//    IO::loadu8_in_01(im_in, std::string(MY_PATH)+argv[1]);

//    im_in.for_all_pixels([&] (ImageGrayd::PixelType &pix)
//    {
//        pix *= pix;
//    });

//    ImageSpectrald modulus, phase;
//    Fourier::fftForwardModulusAndPhase(im_in, modulus, phase);

//    modulus.for_all_pixels([&] (ImageSpectrald::PixelType &pix)
//    {
//        pix = std::abs(pix);
//    });

//    IO::save01_in_u8(modulus, std::string(MY_PATH)+"varSpectrum2_" + name_noext + ".png");

//    return 0;
//}

int test_texton(int argc, char **argv)
{
    if( argc < 3 )
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " <out_path> [source code dependant options]" << std::endl;

        return EXIT_FAILURE;
    }

    QApplication app(argc, argv);
    std::setlocale(LC_ALL,"C");

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

    //warning; this import function turns a texton file into a vizualizable image
    assert(import_texton_from_png(im_texton, MY_PATH+filename_source));

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
    sampler.setNbPoints(300);
    Stamping::StampDiscrete<ImageRGBd> stamp(im_texton);
    stamp.setInterpolationRule(Stamping::StampDiscrete<ImageRGBd>::BILINEAR);

    Stamping::StamperTexton<ImageRGBd> tamponneur(&sampler, &stamp);

    tamponneur.setPeriodicity(false); //if on, the ouput image will be periodic
    tamponneur.setUseMargins(true); //if on, the program will add extra margins to shoot textons on
                                    //when the image is not periodic. It is highly advised you let this on true,
                                    //as the result will lose a lot of energy (variance) otherwise.

    int W=im_in.width(), H=im_in.height();

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

    colored_RPN(im_in, im_rpn, RGB_SPACE, NORMAL, 0, 0, false, true, false, 1.0);

    auto clamp = [&] (ImageRGBd::PixelType &pix) { for(int i=0; i<3; ++i) pix[i] = pix[i] > 1.0 ? 1.0 : (pix[i] < 0.0 ? 0.0 : pix[i]); };
    im_out.for_all_pixels(clamp);
    im_rpn.for_all_pixels(clamp);

    HistogramRGBd histo_in(im_in);
    HistogramRGBd histo_out(im_out);
    HistogramRGBd histo_rpn(im_rpn);

    std::cout << "Mean of input is: (" << std::to_string(histo_in.mean(0)) << ", " << std::to_string(histo_in.mean(1)) << ", " << std::to_string(histo_in.mean(2)) << ")" << std::endl;
    std::cout << "Variance of input is: (" << std::to_string(histo_in.covariance(0, 0)) << ", " << std::to_string(histo_in.covariance(1, 1)) << ", " << std::to_string(histo_in.covariance(2, 2)) << ")" << std::endl;

    std::cout << "==================" << std::endl;

    std::cout << "Mean of output is: (" << std::to_string(histo_out.mean(0)) << ", " << std::to_string(histo_out.mean(1)) << ", " << std::to_string(histo_out.mean(2)) << ")" << std::endl;
    std::cout << "Variance of output is: (" << std::to_string(histo_out.covariance(0, 0)) << ", " << std::to_string(histo_out.covariance(1, 1)) << ", " << std::to_string(histo_out.covariance(2, 2)) << ")" << std::endl;

    std::cout << "==================" << std::endl;

    std::cout << "Mean of rpn is: (" << std::to_string(histo_rpn.mean(0)) << ", " << std::to_string(histo_rpn.mean(1)) << ", " << std::to_string(histo_rpn.mean(2)) << ")" << std::endl;
    std::cout << "Variance of rpn is: (" << std::to_string(histo_rpn.covariance(0, 0)) << ", " << std::to_string(histo_rpn.covariance(1, 1)) << ", " << std::to_string(histo_rpn.covariance(2, 2)) << ")" << std::endl;

    std::cout << "==================" << std::endl;

    std::cout << "Mean of texton is: (" << std::to_string(histo_texton.mean(0)) << ", " << std::to_string(histo_texton.mean(1)) << ", " << std::to_string(histo_texton.mean(2)) << ")" << std::endl;
    std::cout << "Variance of texton is: (" << std::to_string(histo_texton.covariance(0, 0)) << ", " << std::to_string(histo_texton.covariance(1, 1)) << ", " << std::to_string(histo_texton.covariance(2, 2)) << ")" << std::endl;


    double zero[3];
    double one[3];
    for(int i=0; i<3; ++i) {zero[i]=0.0; one[i]=1.0;}

    HistogramRGBBase<int> quantizedHisto;

    quantizedHisto=histo_in.quantize(zero, one, 24);
    quantizedHisto.saveHistogram(std::string(MY_PATH) + "in_" + input_noext + ".csv", 24);

    quantizedHisto=histo_out.quantize(zero, one, 24);
    quantizedHisto.saveHistogram(std::string(MY_PATH) + "out_" + input_noext + "_tn" + ".csv", 24);

    quantizedHisto=histo_rpn.quantize(zero, one, 24);
    quantizedHisto.saveHistogram(std::string(MY_PATH) + "out_" + input_noext + "_rpn" + ".csv", 24);

    if(argc > 3) //add an original texton noise image to compare the results with it
    {
        filename_source=argv[3];
        ImageRGBd im_tn;
        IO::loadu8_in_01(im_tn, MY_PATH + filename_source);

        HistogramRGBd histo_tn(im_tn);

        std::cout << "==================" << std::endl;

        std::cout << "Mean of provided texton noise is: (" << std::to_string(histo_tn.mean(0)) << ", " << std::to_string(histo_tn.mean(1)) << ", " << std::to_string(histo_tn.mean(2)) << ")" << std::endl;
        std::cout << "Variance of provided texton noise is: (" << std::to_string(histo_tn.covariance(0, 0)) << ", " << std::to_string(histo_tn.covariance(1, 1)) << ", " << std::to_string(histo_tn.covariance(2, 2)) << ")" << std::endl;

        quantizedHisto=histo_tn.quantize(zero, one, 24);
        quantizedHisto.saveHistogram(std::string(MY_PATH) + "out_" + input_noext + "_tn_galerne" + ".csv", 24);
    }


    //Saving sample
    //Saving im_out
    IO::save01_in_u8(im_sample, std::string(MY_PATH) + name_noext + "_sample.png");
    IO::save01_in_u8(im_out, std::string(MY_PATH) + name_noext + "_out_tn.png");
    IO::save01_in_u8(im_rpn, std::string(MY_PATH) + input_noext + "_out_rpn.png");
    IO::save01_in_u8(im_texton, std::string(MY_PATH) + name_noext + "_texton.png");

    ImageViewer imgv_in("Source", &app, 0);
    imgv_in.setWindowTitle("Source");
    imgv_in.set_rgb01(im_in.getDataPtr(), im_in.width(),im_in.height(),1);
    imgv_in.show();

    ImageViewer imgv_sample("Sampling", &app, 1);
    imgv_sample.setWindowTitle("Sampling");
    imgv_sample.set_rgb01(im_sample.getDataPtr(), im_sample.width(), im_sample.height(),1);
    imgv_sample.show();

    ImageViewer imgv_rpn("RPN", &app, 2);
    imgv_rpn.setWindowTitle("RPN");
    imgv_rpn.set_rgb01(im_rpn.getDataPtr(), im_rpn.width(),im_rpn.height(),1);
    imgv_rpn.show();

    ImageViewer imgv_out("Texton noise", &app, 3);
    imgv_out.setWindowTitle("Texton noise");
    imgv_out.set_rgb01(im_out.getDataPtr(), im_out.width(), im_out.height(),1);
    imgv_out.show();

    im_texton.for_all_pixels(clamp);

    ImageViewer imgv_texton("Texton", &app, 4);
    imgv_texton.setWindowTitle("Texton");
    imgv_texton.set_rgb01(im_texton.getDataPtr(), im_texton.width(), im_texton.height(),1);
    imgv_texton.show();

    return app.exec();
}

//int test_autoconvolutionSpectrum(int argc, char **argv)
//{
//    if( argc < 3 )
//    {
//        std::cerr << "Usage: " << std::endl;
//        std::cerr << argv[0] << " <out_path> [source code dependant options]" << std::endl;

//        return EXIT_FAILURE;
//    }

//    QApplication app(argc, argv);
//    std::setlocale(LC_ALL,"C");

//    std::string filename_source=argv[1];
//    std::string name_file = IO::remove_path(filename_source);
//    std::string name_noext = IO::remove_ext(name_file);

//    //std::string out_path = argv[1];

//    std::string input_noext=name_noext;

//    filename_source=argv[2];
//    name_file = IO::remove_path(filename_source);
//    name_noext = IO::remove_ext(name_file);

//    //Loading sample
//    //Loading im_in

//    ImageGrayd im_in;
//    IO::loadu8_in_01(im_in, std::string(MY_PATH)+argv[1]);

//    ImageSpectrald modulus, phase;
//    ImageSpectrald mm;
//    modulus.initItk(im_in.width(), im_in.height());
//    mm.initItk(modulus.width(), modulus.height(), true);
//    phase.initItk(im_in.width(), im_in.height());
//    Fourier::fftForwardModulusAndPhase(im_in, modulus, phase);
//    modulus.for_all_pixels([&] (ImageSpectrald::PixelType &p1, int x, int y)
//    {
//        modulus.for_all_pixels([&] (ImageSpectrald::PixelType &p2)
//        {
//            mm.pixelAbsolute(x, y) += p1 * p2;
//        });
//    });

//    ImageViewer imgv_in("Source", &app, 0);
//    imgv_in.set_gray01(im_in.getDataPtr(), im_in.width(), im_in.height(), 1);
//    imgv_in.show();

//    ImageViewer imgv_modulus("Input PSD", &app, 1);
//    imgv_modulus.set_gray01(modulus.getDataPtr(), modulus.width(), modulus.height(), 1);
//    imgv_modulus.show();

//    ImageViewer imgv_out("Output PSD", &app, 2);
//    imgv_out.set_gray01(mm.getDataPtr(), mm.width(), mm.height(), 1);
//    imgv_out.show();

//    return app.exec();
//}

int main( int argc, char **argv )
{
    //test_genet(argc, argv);
    test_texton(argc, argv);
    //test_autoconvolutionSpectrum(argc, argv);
}
