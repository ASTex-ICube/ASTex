#include <ASTex/easy_io.h>
#include "pnf.h"

using namespace ASTex;

int main()
{
    Noise<T> noise(100,10.,20.);
    Color_map<T> cm;

//  palette a
    cm.add_color(0,Color(1,1,0));
    cm.add_color(40,Color(1,0,0));
    cm.add_color(59,Color(0,0,0));
    cm.add_color(60,Color(1,1,1));
    cm.add_color(100,Color(1,1,1));

//  palette b
//    cm.add_color(0, Color(0,0,0));
//    cm.add_color(4, Color(0,0,0));
//    cm.add_color(5, Color(1,1,1));
//    cm.add_color(7, Color(1,1,1));
//    cm.add_color(9, Color(199./255., 139./255., 105./255.));
//    cm.add_color(10, Color(199./255., 139./255., 105./255.));

//  palette c
//    cm.add_color(0, Color(0., 0., 1.));
//    cm.add_color(1, Color(1., 0., 0.));

//  palette d
//    cm.add_color(0,Color(0,0,1));
//    cm.add_color(1,Color(0,1,0));
//    cm.add_color(2,Color(1,0,0));

    cm.export_courbe(TEMPO_PATH + "data.txt");
//    cm.export_img_palette(512, TEMPO_PATH + "palette.png");


    //filtrage color map
    auto start_chrono = std::chrono::system_clock::now();

    cm.filter(512,512,200,0.3);

    std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "filtrage de la color map timing: " << elapsed_seconds.count() << " s." << std::endl;

    ImageRGB<T> c0_ = cm.get_filtered();
    IO::save01_in_u8(c0_,TEMPO_PATH + "color_map_filtered.png");

//    ImageRGB<T> c0_;
//    IO::loadu8_in_01(c0_,TEMPO_PATH+ "color_map_filtered.png");
//    cm.set_filtered(c0_,T(0.5));

    Vec2 w_size(2*256,2*256);
    Vec2 im_size(512,512);

    //unfiltered noise

    start_chrono = std::chrono::system_clock::now();
    ImageGray<T> noise_unfiltered = compute_noise_unfiltered(w_size,im_size,noise);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe unfiltering timing: " << elapsed_seconds.count() << " s." << std::endl;

    IO::save01_in_u8(noise_unfiltered,TEMPO_PATH + "noise_unfiltered_sum_cosines_2x2.png");

    //filtered noise

    start_chrono = std::chrono::system_clock::now();
    ImageGray<T> noise_filtered = compute_noise_filtered(w_size,im_size,noise, 10);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe noise filtering timing: " << elapsed_seconds.count() << " s." << std::endl;

    IO::save01_in_u8(noise_unfiltered,TEMPO_PATH + "noise_filtered_sum_cosines_2x2.png");

    //unfiltered noise mapped
    start_chrono = std::chrono::system_clock::now();
    ImageRGB<T> unfiltered = compute_unfiltered_IMG(w_size,im_size,noise,cm);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe unfiltering mapped timing: " << elapsed_seconds.count() << " s." << std::endl;

    IO::save01_in_u8(unfiltered,TEMPO_PATH + "noise_mapped_unfiltered_sum_cosines_2x2.png");

    //ground truth

    start_chrono = std::chrono::system_clock::now();
    ImageRGB<T> ground_truth = compute_ground_truth_IMG(w_size,im_size,noise,cm,10);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe ground_truth timing: " << elapsed_seconds.count() << " s." << std::endl;

    IO::save01_in_u8(ground_truth,TEMPO_PATH + "noise_mapped_ground_truth_sum_cosines_2x2.png");

    //naive filtering
    start_chrono = std::chrono::system_clock::now();
    ImageRGB<T> naive = compute_naive_filter_IMG(w_size,im_size,noise,cm,10);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe naive filtering timing: " << elapsed_seconds.count() << " s." << std::endl;

    IO::save01_in_u8(naive,TEMPO_PATH + "noise_mapped_naive_filtered_sum_cosines_2x2.png");

    //good filtering
    start_chrono = std::chrono::system_clock::now();
    ImageRGB<T> filtered = compute_good_filter_IMG(w_size,im_size,noise,cm,10);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe good filtering timing: " << elapsed_seconds.count() << " s." << std::endl;

    IO::save01_in_u8(filtered,TEMPO_PATH + "noise_mapped_good_filtered_sum_cosines_2x2.png");

    return 0;
}
