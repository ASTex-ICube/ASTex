//
// Created by grenier on 26/06/23.
//
#include <ASTex/easy_io.h>
#include <ASTex/image_rgb.h>
#include "ASTex/Noises/Gabor.h"


int main(){

// ---------------------------------------------------------------------------
// load d'image
    int cm_size = 256;

    ImageRGBu8 cm_(cm_size, cm_size);
//    cm_.load("/home/grenier/Documents/ASTex_fork/results/T_analysis/cm_in.png");
    cm_.load("/home/grenier/Documents/ASTex_fork/results/T_analysis/cm_out.png");

// ---------------------------------------------------------------------------
// génération noise

    int resolution = 256;
    int img_size = 512;

    float F_0_ = 0.04;
    float omega_0_ = 0.;//M_PI;

    float number_of_impulses_per_kernel = 64.0;
    unsigned period = 128; // non utilisé

    float K_ = 1.0;
    float a_ = 0.02;

    unsigned random_offset_ = 954248632;

    noise noise_1(K_,
                 a_,
                 F_0_,
                 omega_0_,
                 number_of_impulses_per_kernel,
                 period,
                 random_offset_,
                 1);

    noise noise_2(K_,
                 a_,
                 F_0_,
                 omega_0_,
                 number_of_impulses_per_kernel,
                 period,
                 random_offset_,
                 2);

    ImageGrayu8 noise1_ = storing_noise(resolution, img_size, noise_1);
    ImageGrayu8 noise2_ = storing_noise(resolution, img_size, noise_2);

    noise1_.save("/home/grenier/Documents/ASTex_fork/results/T_analysis/gabor1.png");
    noise2_.save("/home/grenier/Documents/ASTex_fork/results/T_analysis/gabor2.png");



// ---------------------------------------------------------------------------
// composition
    ImageRGBu8 result(img_size, img_size);
    result.parallel_for_all_pixels([&] (typename ImageRGBu8::PixelType& P, int x, int y)
                                   {
                                       double n1 = noise1_.pixelAbsolute(x,y);
                                       double n2 = noise2_.pixelAbsolute(x,y);

                                       int id_n1 = std::clamp(int(n1), 0, img_size);
                                       int id_n2 = std::clamp(int(n2), 0, img_size);

                                       P = ImageRGBu8::PixelType(cm_.pixelAbsolute(id_n1, id_n2));
                                   });


// ---------------------------------------------------------------------------
// result
//    result.save("/home/grenier/Documents/ASTex_fork/results/T_analysis/composition_in.png");
    result.save("/home/grenier/Documents/ASTex_fork/results/T_analysis/composition_out.png");

    return 0;
}

