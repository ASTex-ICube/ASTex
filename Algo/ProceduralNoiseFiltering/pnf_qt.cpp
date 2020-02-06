#include "texture_noise.h"
#include "pnf.h"
#include <ASTex/color_map.h>
#include <ASTex/rpn_utils.h>
#include "imageviewer.h"

using namespace ASTex;

int main(int argc, char **argv)
{
    QApplication app(argc, argv);


    ImageSpectrald psd;
    IO::EXR::load(psd, TEMPO_PATH + "spectra/donut_black.exr");

//    ImageGray<T> example_noise;
//    IO::loadu8_in_01(example_noise, TEMPO_PATH + "gray_png/gc04.png");

//    ImageSpectrald psd, phase;
//    Fourier::fftForwardModulusAndPhase(example_noise, psd, phase);

    TextureNoise<T> texture_noise(psd);
    Color_map<T> cm;

    ImageRGB<T> c0_;
    IO::loadu8_in_01(c0_,TEMPO_PATH+ "color_map_filtered.png");
    cm.set_filtered(c0_,0.5);

    Vec2 w_size(64,64);
    Vec2 im_size(512,512);


    // compute

    ImageGray<T> noise = texture_noise.getNoise();


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


    //show window

//    ImageViewer iv_en("Example noise");
//    iv_en.update(example_noise);
//    iv_en.show();

    ImageViewer iv_psd("PSD");
    iv_psd.update(psd);
    iv_psd.show();

    ImageViewer iv_noise("Noise");
    iv_noise.update(noise);
    iv_noise.show();

    IO::save01_in_u8(noise, TEMPO_PATH + "noise_PSD.png");

    ImageViewer iv_cm("Color map");
    iv_cm.update(cm.get_filtered());
    iv_cm.show();

    ImageViewer iv_noise_cm("Noise color mapped");
    iv_noise_cm.update(im_noise_cm);
    iv_noise_cm.show();

    ImageViewer iv_ground_truth("Ground truth");
    iv_ground_truth.update(im_ground_truth);
    iv_ground_truth.show();

    ImageViewer iv_noise_cm_naive_filter("Noise color mapped with naive filtering");
    iv_noise_cm_naive_filter.update(im_noise_cm_naive_filter);
    iv_noise_cm_naive_filter.show();

    ImageViewer iv_noise_cm_good_filter("Noise color mapped with good filtering");
    iv_noise_cm_good_filter.update(im_noise_cm_good_filter);
    iv_noise_cm_good_filter.show();

    return app.exec();
}
