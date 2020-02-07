#include <ASTex/easy_io.h>
#include "pnf.h"

using namespace ASTex;

int main()
{
    Noise<T> noise(1000,10.,20.);
    Color_map<T> cm;

//  palette 1
    cm.add_color(0,Color(1,1,0));
    cm.add_color(40,Color(1,0,0));
    cm.add_color(59,Color(0,0,0));
    cm.add_color(60,Color(1,1,1));
    cm.add_color(100,Color(1,1,1));


//  palette 2
//    cm.add_color(0,Color(0,0,1));
//    cm.add_color(1,Color(0,1,0));
//    cm.add_color(2,Color(1,0,0));

//  palette 3
//    cm.add_color(0, Color(154./255., 173./255., 213./255.));
//    cm.add_color(1, Color(196./255., 192./255., 144./255.));

//  palette 4
//    cm.add_color(0, Color(0,0,0));
//    cm.add_color(4, Color(0,0,0));
//    cm.add_color(5, Color(1,1,1));
//    cm.add_color(7, Color(1,1,1));
//    cm.add_color(9, Color(199./255., 139./255., 105./255.));
//    cm.add_color(10, Color(199./255., 139./255., 105./255.));

    cm.export_courbe(TEMPO_PATH + "data.txt");
    cm.export_img_palette(512, TEMPO_PATH + "palette.png");


    //filtrage color map
    auto start_chrono = std::chrono::system_clock::now();

    cm.filter(512,512,200,0.5);

    std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "filtrage de la color map timing: " << elapsed_seconds.count() << " s." << std::endl;

    ImageRGB<T> c0_ = cm.get_filtered();
    IO::save01_in_u8(c0_,TEMPO_PATH + "color_map_filtered.png");

//    ImageRGB<T> c0_;
//    IO::loadu8_in_01(c0_,TEMPO_PATH+ "color_map_filtered.png");
//    cm.set_filtered(c0_,T(0.5));

    Vec2 w_size(32,32);
    Vec2 im_size(512,512);

    //unfilterd noise
    start_chrono = std::chrono::system_clock::now();
    ImageRGB<T> unfiltered = compute_unfiltered_IMG(w_size,im_size,noise,cm);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe unfiltering timing: " << elapsed_seconds.count() << " s." << std::endl;

    IO::save01_in_u8(unfiltered,TEMPO_PATH + "noise_unfiltered_sum_cosines_32x32.png");

    //ground truth

    start_chrono = std::chrono::system_clock::now();
    ImageRGB<T> ground_truth = compute_ground_truth_IMG(w_size,im_size,noise,cm,10);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe ground_truth timing: " << elapsed_seconds.count() << " s." << std::endl;

    IO::save01_in_u8(ground_truth,TEMPO_PATH + "ground_truth_sum_cosines_32x32.png");

    //naive filtering
    start_chrono = std::chrono::system_clock::now();
    ImageRGB<T> naive = compute_naive_filter_IMG(w_size,im_size,noise,cm,10);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe naive filtering timing: " << elapsed_seconds.count() << " s." << std::endl;

    IO::save01_in_u8(naive,TEMPO_PATH + "noise_naive_filtered_sum_cosines_32x32.png");

    //good filtering
    start_chrono = std::chrono::system_clock::now();
    ImageRGB<T> filtered = compute_good_filter_IMG(w_size,im_size,noise,cm,10);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe good filtering timing: " << elapsed_seconds.count() << " s." << std::endl;

    IO::save01_in_u8(filtered,TEMPO_PATH + "noise_filtered_sum_cosines_32x32.png");

    return 0;
}
