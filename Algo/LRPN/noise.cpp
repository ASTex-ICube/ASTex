
#include "io.h"
#include "fourier.h"
#include "itkJPEGImageIO.h"
#include "mode_seeking.h"
#include "output_geoffrey.h"
#include "image_treatment.h"
#include "utilities.h"
#include "tests.h"
#include "colorspace_filters.h"


using namespace ASTex;


/**
 * @brief RPnoise_from_spectrum random phase noise
 * @param inputfile a spectrum (encoded as grayscale image)
 * @param outputfile is a patchwork image containing spectrum + output (pure random phase)
 * @return error code
 */
int RPnoise_from_spectrum(std::string inputfile, std::string outputfile, int output_w, int output_h, int blending_size)
{
    // COLLECTION OF IMAGES
    image_collector collec;

    // LOAD MODULUS
    ImageGrayd modulus;
    load_spectrum(modulus,inputfile);
    const int im_w = modulus.width();
    const int im_h = modulus.height();

    // MODIFY INPUT FOR PRINT
    ImageGrayd input;
    input.initItk(im_w,im_h,true);
    input.setCenter(im_w/2,im_h/2);

//    std::cout << distance_spectrum_to_spectrum(modulus, modulus) << std::endl;
//    std::cout << distance_spectrum_to_spectrum(modulus, input) << std::endl;

    for (int x=0; x< im_w; ++x)
    {
        for (int y=0; y<im_h; ++y)
        {
            const double v = modulus.pixelAbsolute(x,y);
            input.pixelAbsolute(x,y) =  1 - std::exp(- v*2.0);
        }
    }
    collec.add(input); // add spectrum to collection
//    std::cout << distance_spectrum_to_spectrum(modulus, input) << std::endl;

    // RESULT
    ImageGrayd result;
    result.initItk(output_w,output_h,true);
    RPnoise_mosaic(modulus,result,blending_size);

    collec.add(result,0.0); // add result to collection
    save(collec.collect(),outputfile, 0.0);

    return EXIT_SUCCESS;
}


/**
 * @brief main function
 * @return error code
 */
int main()
{
//    return RPnoise_from_spectrum(TEMPO_PATH+"spectrum.png",TEMPO_PATH+"/result.png",512,512,4);
    return testSpectrumExtraction3(TEMPO_PATH+"noise.png",TEMPO_PATH+"/mask_0.png",TEMPO_PATH+"/mask_1.png",TEMPO_PATH+"/result.png");
//    return testSpectrumExtraction2(TEMPO_PATH+"refnoise.png",TEMPO_PATH+"/result.png");
//    return RPnoise(TEMPO_PATH+"refnoise.png",TEMPO_PATH+"/result.png");
//  return testMask(TEMPO_PATH+"refnoise.png",TEMPO_PATH+"/result.png");
//  return testPhase(TEMPO_PATH+"refnoise.png","/tmp/result.png");
//	return testAnimation(TEMPO_PATH+"refnoise.png","/tmp/result_");
}
