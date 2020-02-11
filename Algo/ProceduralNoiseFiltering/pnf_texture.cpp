#include "texture_noise.h"
#include "pnf.h"
#include <ASTex/color_map.h>
#include <ASTex/rpn_utils.h>

using namespace ASTex;

int main()
{
    ImageGray<T> noise;
    IO::loadu16_in_01(noise, TEMPO_PATH + "noise/voronoi_repeat_non_gauss.png");

//    ImageSpectrald psd;
//    IO::loadu8_in_01(psd, TEMPO_PATH + "spectra/donut_black.png");

//    ImageGray<T> example_noise;
//    IO::loadu8_in_01(example_noise, TEMPO_PATH + "gray_png/gc12.png");

//    ImageSpectrald psd, phase;
//    Fourier::fftForwardModulusAndPhase(example_noise, psd, phase);

    TextureNoise<T> texture_noise;
    texture_noise.setNoise(noise);
    Color_map<T> cm;

    ImageRGB<T> c0_;
    IO::loadu8_in_01(c0_,TEMPO_PATH+ "color_map_filtered.png");
    cm.set_filtered(c0_,0.3);

    Vec2 w_size(noise.width()*0.125, noise.height()*0.125);
    Vec2 im_size(1024,1024);


    // compute

//    ImageGray<T> noise = texture_noise.getNoise();


    // noise unfiltered
    auto start_chrono = std::chrono::system_clock::now();

    ImageRGB<T> im_noise_cm = compute_unfiltered_IMG(w_size, im_size, texture_noise, cm);

    std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe unfiltering timing: " << elapsed_seconds.count() << " s." << std::endl;


    // ground truth
    start_chrono = std::chrono::system_clock::now();

    ImageRGB<T> im_ground_truth = compute_ground_truth_IMG(w_size, im_size, texture_noise, cm);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe ground_truth timing: " << elapsed_seconds.count() << " s." << std::endl;


    // naive filter
    start_chrono = std::chrono::system_clock::now();

    ImageRGB<T> im_noise_cm_naive_filter = compute_naive_filter_IMG(w_size, im_size, texture_noise, cm);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe naive filtering timing: " << elapsed_seconds.count() << " s." << std::endl;

    //good filter
    start_chrono = std::chrono::system_clock::now();

    ImageRGB<T> im_noise_cm_good_filter = compute_good_filter_IMG(w_size, im_size, texture_noise, cm);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe good filtering timing: " << elapsed_seconds.count() << " s." << std::endl;

//    IO::save01_in_u8(noise, TEMPO_PATH + "texture_noise.png");
//    IO::save01_in_u8(im_noise_cm,TEMPO_PATH + "texture_noise_unfilered.png");
//    IO::save01_in_u8(im_ground_truth,TEMPO_PATH + "texture_ground_truth.png");
//    IO::save01_in_u8(im_noise_cm_naive_filter,TEMPO_PATH + "texture_noise_naive_filtering.png");
//    IO::save01_in_u8(im_noise_cm_good_filter,TEMPO_PATH + "texture_noise_goood_filtering.png");

//    IO::save01_in_u16(noise, TEMPO_PATH + "texture2_noise.png");
    IO::save01_in_u16(im_noise_cm,TEMPO_PATH + "texture2_noise_unfilered.png");
    IO::save01_in_u16(im_ground_truth,TEMPO_PATH + "texture2_ground_truth.png");
    IO::save01_in_u16(im_noise_cm_naive_filter,TEMPO_PATH + "texture2_noise_naive_filtering.png");
    IO::save01_in_u16(im_noise_cm_good_filter,TEMPO_PATH + "texture2_noise_goood_filtering.png");

    return 0;
}
