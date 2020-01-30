#include <ASTex/easy_io.h>
#include "pnf.h"

using namespace ASTex;

int main()
{
    Noise<T> noise;
    Color_map<T> cm;
    cm.add_color(0,Color(1,1,0));
    cm.add_color(40,Color(1,0,0));
    cm.add_color(59,Color(0,0,0));
    cm.add_color(60,Color(1,1,1));
    cm.add_color(100,Color(1,1,1));


    //cm.export_palette(TEMPO_PATH + "palette.gnu");
    cm.export_courbe(TEMPO_PATH + "data.txt");


    //filtrage color map
    auto start_chrono = std::chrono::system_clock::now();

    cm.filter(512,512,200,T(1)/T(6));

    std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "filtrage de la color map timing: " << elapsed_seconds.count() << " s." << std::endl;

    ImageRGB<T> c0_ = cm.get_filtered();
    IO::save01_in_u8(c0_,TEMPO_PATH + "color_map_filtered.png");

//    ImageRGB<T> c0_;
//    IO::loadu8_in_01(c0_,TEMPO_PATH+ "color_map_filtered.png");
//    cm.set_filtered(c0_,T(0.5));

    Vec2 w_size(512,512);
    Vec2 im_size(512,512);

    //unfilterd noise
    start_chrono = std::chrono::system_clock::now();
    ImageRGB<T> unfiltered = compute_unfiltered_IMG(w_size,im_size,noise,cm);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe unfiltering timing: " << elapsed_seconds.count() << " s." << std::endl;

    IO::save01_in_u8(unfiltered,TEMPO_PATH + "noise_unfiltered.png");

    //naive filtering
    start_chrono = std::chrono::system_clock::now();
    ImageRGB<T> naive = compute_naive_filter_IMG(w_size,im_size,noise,cm,10);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe naive filtering timing: " << elapsed_seconds.count() << " s." << std::endl;

    IO::save01_in_u8(naive,TEMPO_PATH + "noise_naive_filtered.png");

    //good filtering
    start_chrono = std::chrono::system_clock::now();
    ImageRGB<T> filtered = compute_good_filter_IMG(w_size,im_size,noise,cm,10);

    elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe good filtering timing: " << elapsed_seconds.count() << " s." << std::endl;

    IO::save01_in_u8(filtered,TEMPO_PATH + "noise_filtered.png");

    return 0;
}
