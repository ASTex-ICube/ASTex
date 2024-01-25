//
// Created by grenier on 26/06/23.
//


#include <cmath>
#include <ASTex/easy_io.h>
#include <ASTex/image_rgb.h>
#include "ASTex/histogram.h"

#include "ASTex/fourier.h"




void fourier_analysis(ImageGrayd input, std::string name){
    ImageSpectrald module;
    ImageSpectrald phase;
    ImageGrayd result;

    Fourier::fftForwardModulusAndPhase(input, module, phase, false);
    Fourier::fftInverseModulusAndPhase(module, phase, result, false);

//    IO::save01_in_u8(input,"/home/grenier/Documents/ASTex_fork/results/fourier/res/" + name + "_input_noise.png");
    IO::save_spectrum(module, "/home/grenier/Documents/ASTex_fork/results/fourier/res/" + name + "_modulus.png");
    IO::save_phase(phase, "/home/grenier/Documents/ASTex_fork/results/fourier/res/" + name + "_phase.png");
    IO::save01_in_u8(result,"/home/grenier/Documents/ASTex_fork/results/fourier/res/" + name + "_inv_noise.png");
}






int main(){

// ---------------------------------------------------------------------------
// load d'une image

    ImageGrayd image_(1024, 1024);
    IO::loadu8_in_01(image_, "/home/grenier/Documents/ASTex_fork/results/fourier/cm_complexe_2.png");

// ---------------------------------------------------------------------------
//// test
//    image_.parallel_for_all_pixels([&] (typename ImageGrayd::PixelType& P, int x, int y) // cm
//                                     {
//                                         double lvl = 0.12;
//                                         double pix = image_.pixelAbsolute(x,y);
//
//                                         if(pix < lvl){
//                                             P = ImageGrayd::PixelType(0.);
//                                         }
//                                         else{
//                                             P = ImageGrayd::PixelType(1.);
//                                         }
//                                     });


// ---------------------------------------------------------------------------
// fourier

    fourier_analysis(image_, "res_t");

    return 0;
}


