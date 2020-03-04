#include "texture_noise.h"
#include "pnf.h"
#include "color_map.h"
#include "gaussian_transfer.h"
#include "histogram.h"
#include <ASTex/rpn_utils.h>

using namespace ASTex;

int main()
{
//    using IMG = ImageRGBf;
//    IMG input;
//    IO::loadu8_in_01(input, TEMPO_PATH + "input.png");

//    int nb_bins = 256;
//    float ymax = 0.025f;
//    Histogram<IMG> h;
//    h.computeHisto(input,nb_bins);
//    h.exportHisto(TEMPO_PATH+"histo_input", ymax);

//    ImageRGBf output;
//    IO::loadu8_in_01(output, TEMPO_PATH + "output.png");

//    IMG input_t(input.width(), input.height());
//    Gaussian_transfer::ComputeTinput(input,input_t);

//    h.computeHisto(input_t,nb_bins);
//    h.exportHisto(TEMPO_PATH+"histo_input_t", ymax);

//    IO::save01_in_u8(input_t,TEMPO_PATH + "input_t.png");

//    IMG lut(128,1);
//    Gaussian_transfer::ComputeinvT(input,lut);

//    IO::save01_in_u8(lut,TEMPO_PATH + "lut.png");

//    IMG input_t_t_1(input_t.width(), input_t.height());
//    input_t_t_1.parallel_for_all_pixels([&](IMG::PixelType &pix,int x,int y){
//       auto p = input_t.pixelAbsolute(x,y);
//       //auto p = output.pixelAbsolute(x,y);
//       float r = p.GetRed() * (lut.width()-1);
//       float g = p.GetGreen() * (lut.width()-1);
//       float b = p.GetBlue() * (lut.width()-1);

//       r = lut.pixelAbsolute((int)r,0).GetRed();
//       g = lut.pixelAbsolute((int)g,0).GetGreen();
//       b = lut.pixelAbsolute((int)b,0).GetBlue();

//       pix =  IMG::itkPixel(r, g, b);

////       pix = lut.pixelAbsolute(int(p * (lut.width()-1)),0);
//    });

//    h.computeHisto(input_t_t_1,nb_bins);
//    h.exportHisto(TEMPO_PATH+"histo_input_t_t_1", ymax);

//    IO::save01_in_u16(input_t_t_1, TEMPO_PATH + "input_t_t_1.png");

//    return EXIT_SUCCESS;

    ImageGray<T> noise;
    IO::loadu16_in_01(noise, TEMPO_PATH + "noise/voronoi_repeat_non_gauss.png");
    ImageGray<T> noise_gauss(noise.width(), noise.height());
	Gaussian_transfer<ImageGray<T>>::ComputeTinput(noise,noise_gauss);

//    ImageSpectrald psd;
//    IO::EXR::load(psd, TEMPO_PATH + "spectra/donut_black.exr");

//    ImageGray<T> example_noise;
//    IO::loadu8_in_01(example_noise, TEMPO_PATH + "gray_png/gc12.png");

//    ImageSpectrald psd, phase;
//    Fourier::fftForwardModulusAndPhase(example_noise, psd, phase);

    TextureNoise<T> texture_noise;
    texture_noise.setNoise(noise);
    Color_map<T> cm;

    ImageRGB<T> c0_;
    IO::loadu8_in_01(c0_,TEMPO_PATH+ "color_map_filtered1.png");
    cm.set_filtered(c0_,0.5f);

    Vec2 w_size(noise.width()*16, noise.height()*16);
    Vec2 im_size(noise.width(),noise.height());


    // compute

    //ImageGray<T> noise = texture_noise.getNoise();

    //noise unfilered
    auto start_chrono = std::chrono::system_clock::now();

    ImageGray<T> im_noise_unfiltered = compute_noise_unfiltered(w_size, im_size, texture_noise);

    std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe noise unfiltered timing: " << elapsed_seconds.count() << " s." << std::endl;

    //noise filtered
    start_chrono = std::chrono::system_clock::now();

    ImageGray<T> im_noise_filtered = compute_noise_filtered(w_size, im_size, texture_noise);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe noise filtered timing: " << elapsed_seconds.count() << " s." << std::endl;


    // noise mapped unfiltered
    start_chrono = std::chrono::system_clock::now();

    ImageRGB<T> im_noise_cm = compute_unfiltered_IMG(w_size, im_size, texture_noise, cm);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe noise mapped unfiltered timing: " << elapsed_seconds.count() << " s." << std::endl;


    // ground truth
    start_chrono = std::chrono::system_clock::now();

    ImageRGB<T> im_ground_truth = compute_ground_truth_IMG(w_size, im_size, texture_noise, cm);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe ground truth timing: " << elapsed_seconds.count() << " s." << std::endl;


    // naive filter
    start_chrono = std::chrono::system_clock::now();

    ImageRGB<T> im_noise_cm_naive_filter = compute_naive_filter_IMG(w_size, im_size, texture_noise, cm);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe noise mapped naive filtering timing: " << elapsed_seconds.count() << " s." << std::endl;

    //good filter
    start_chrono = std::chrono::system_clock::now();

    ImageRGB<T> im_noise_cm_good_filter = compute_good_filter_IMG(w_size, im_size, texture_noise, cm);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe noise mapped good filtering timing: " << elapsed_seconds.count() << " s." << std::endl;

//    IO::save01_in_u8(noise, TEMPO_PATH + "texture_noise_example.png");
//    IO::save01_in_u8(noise_gauss, TEMPO_PATH + "texture_noise_gauss_example.png");
//    IO::save01_in_u8(im_noise_filtered,TEMPO_PATH + "texture_noise_scalar_filered.png");
//    IO::save01_in_u8(im_noise_unfiltered,TEMPO_PATH + "texture_noise_scalar_unfilered.png");
//    IO::save01_in_u8(im_noise_cm,TEMPO_PATH + "texture_noise_mapped_unfilered.png");
//    IO::save01_in_u8(im_ground_truth,TEMPO_PATH + "texture_noise_mapped_ground_truth.png");
//    IO::save01_in_u8(im_noise_cm_naive_filter,TEMPO_PATH + "texture_noise_mapped_naive_filtering.png");
//    IO::save01_in_u8(im_noise_cm_good_filter,TEMPO_PATH + "texture_noise_mapped_goood_filtering.png");

    IO::save01_in_u16(noise, TEMPO_PATH + "texture_noise_example.png");
    IO::save01_in_u16(noise_gauss, TEMPO_PATH + "texture_noise_gauss_example.png");
    IO::save01_in_u16(im_noise_unfiltered,TEMPO_PATH + "texture_noise_scalar_unfilered.png");
    IO::save01_in_u16(im_noise_filtered,TEMPO_PATH + "texture_noise_scalar_filered.png");
    IO::save01_in_u16(im_noise_cm,TEMPO_PATH + "texture_noise_mapped_unfilered.png");
    IO::save01_in_u16(im_ground_truth,TEMPO_PATH + "texture_noise_mapped_ground_truth.png");
    IO::save01_in_u16(im_noise_cm_naive_filter,TEMPO_PATH + "texture_noise_mapped_naive_filtering.png");
    IO::save01_in_u16(im_noise_cm_good_filter,TEMPO_PATH + "texture_noise_mapped_goood_filtering.png");

    //Histogram<ImageGray<T>> gray;
	Histogram2<ImageRGB<T>> rgb;

    rgb.computeHisto(im_ground_truth,256);
    rgb.exportHisto(TEMPO_PATH + "ground_truth", 1.0);

    rgb.computeHisto(im_noise_cm_good_filter,256);
    rgb.exportHisto(TEMPO_PATH + "good_filter", 1.0);

    return EXIT_SUCCESS;
}
