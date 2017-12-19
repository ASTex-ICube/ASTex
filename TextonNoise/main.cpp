#include <stdlib.h>
#include "ASTex/easy_io.h"
#include "ASTex/utils.h"
#include "texton_cpu.h"
#include "histogram.h"
#include "rpn_utils.h"

#define MY_PATH "/home/nlutz/img/texton/"

using namespace ASTex;

int main( int argc, char **argv )
{
    if( argc < 3 )
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " <out_path> [source code dependant options]" << std::endl;

        return EXIT_FAILURE;
    }

//    ImageRGBd im_1, im_2, im_3, im_4;
//    IO::loadu8_in_01(im_1, std::string(MY_PATH)+argv[1]);
//    IO::loadu8_in_01(im_2, std::string(MY_PATH)+argv[2]);
//    IO::loadu8_in_01(im_3, std::string(MY_PATH)+argv[3]);
//    IO::loadu8_in_01(im_4, std::string(MY_PATH)+argv[4]);

//    HistogramRGBd histo_1(im_1);
//    HistogramRGBd histo_2(im_2);
//    HistogramRGBd histo_3(im_3);
//    HistogramRGBd histo_4(im_4);

//    double zero[3];
//    double one[3];
//    for(int i=0; i<3; ++i) {zero[i]=0.0; one[i]=1.0;}

//    HistogramRGBBase<int> quantizedHisto;

//    quantizedHisto=histo_1.quantize(zero, one, 32);
//    quantizedHisto.saveHistogram(std::string(MY_PATH) + "texton1.csv", 32);

//    quantizedHisto=histo_2.quantize(zero, one, 32);
//    quantizedHisto.saveHistogram(std::string(MY_PATH) + "texton2.csv", 32);

//    quantizedHisto=histo_3.quantize(zero, one, 32);
//    quantizedHisto.saveHistogram(std::string(MY_PATH) + "out1.csv", 32);

//    quantizedHisto=histo_4.quantize(zero, one, 32);
//    quantizedHisto.saveHistogram(std::string(MY_PATH) + "out2.csv", 32);

    std::string filename_source=argv[1];
    std::string name_file = IO::remove_path(filename_source);
    std::string name_noext = IO::remove_ext(name_file);

    //std::string out_path = argv[1];

    for(int i=2; i<argc; ++i)
    {
        std::string input_noext=name_noext;

        filename_source=argv[i];
        name_file = IO::remove_path(filename_source);
        name_noext = IO::remove_ext(name_file);

        //Loading sample
        //Loading im_in

        ImageRGBd im_in, im_out, im_sample, im_texton;
        IO::loadu8_in_01(im_in, std::string(MY_PATH)+argv[1]);
        IO::loadu8_in_01(im_texton, MY_PATH+filename_source);

        //Testing

        RandomSampling sampler;
        std::vector<vec2> pointArray = sampler.Generate(2);


        TextonStamper tamponneur(pointArray, im_texton);

        tamponneur.setPeriodicity(false);
        tamponneur.setUseMargins(true);

        int W=im_in.width(), H=im_in.height();

        im_out = tamponneur.generate(W, H);

        im_sample.initItk(W, H, true);

        for(std::vector<vec2>::iterator it=pointArray.begin(); it!=pointArray.end(); ++it)
        {
            int i = im_sample.width() * (*it)[0]; //i & j: single point coordinates in im_out
            int j = im_sample.height() * (*it)[1];

            im_sample.pixelAbsolute(i, j)=1.0;
        }

        HistogramRGBd histo_in(im_in);
        HistogramRGBd histo_out(im_out);

        std::cout << "Variance of input is: (" << std::to_string(histo_in.covariance(0, 0)) << ", " << std::to_string(histo_in.covariance(1, 1)) << ", " << std::to_string(histo_in.covariance(2, 2)) << ")" << std::endl;
        std::cout << "Variance of output is: (" << std::to_string(histo_out.covariance(0, 0)) << ", " << std::to_string(histo_out.covariance(1, 1)) << ", " << std::to_string(histo_out.covariance(2, 2)) << ")" << std::endl;

        double zero[3];
        double one[3];
        for(int i=0; i<3; ++i) {zero[i]=0.0; one[i]=1.0;}

        HistogramRGBBase<int> quantizedHisto;

        quantizedHisto=histo_in.quantize(zero, one, 24);
        quantizedHisto.saveHistogram(std::string(MY_PATH) + "in.csv", 24);

        quantizedHisto=histo_out.quantize(zero, one, 24);
        quantizedHisto.saveHistogram(std::string(MY_PATH) + "out_" + input_noext + "_noperio" + ".csv", 24);

        im_out.for_all_pixels([&] (ImageRGBd::PixelType &pix) {
            for(int i=0; i<3; ++i)
            {
                pix[i] = pix[i] > 1.0 ? 1.0 : (pix[i] < 0.0 ? 0.0 : pix[i]);
            }
        });

        //Saving sample
        //Saving im_out

        //IO::save01_in_u8(im_out, MY_PATH+out_path + name_noext + ".png");
        //IO::save01_in_u8(im_in, std::string(MY_PATH) + "test_input.png");
        IO::save01_in_u8(im_out, std::string(MY_PATH) + "test_textonNoise.png");
        IO::save01_in_u8(im_sample, std::string(MY_PATH) + "test_sample.png");
        IO::save01_in_u8(im_texton, std::string(MY_PATH) + "test_texton.png");
    }




    return EXIT_SUCCESS;
}

