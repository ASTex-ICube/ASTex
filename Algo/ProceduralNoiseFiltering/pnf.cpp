#include <ASTex/noise.h>
#include <ASTex/easy_io.h>
#include <ASTex/color_map.h>
#include <Eigen/Eigen>
#include <ASTex/mipmap.h>
#include <ASTex/rpn_utils.h>

using namespace ASTex;
using T = double;
using Vec2 = Eigen::Matrix<T,2,1>;
using Mat22 = Eigen::Matrix<T,2,2>;
using Color = Color_map<T>::Color;


inline ImageRGB<T> computeNoiseIMG(const Vec2 & w_size, const Vec2 &im_size,const Noise<T> &n,const Color_map<T> &cm){
    Mat22 borns;
    borns << w_size(0), 0, 0, w_size(1);
    Vec2 center = Vec2(borns(0),borns(3)) / T(2);

    //footprint pixel
    Vec2 s(w_size(0) / im_size(0), w_size(1) / im_size(1));

    std::random_device rd;
    std::mt19937 gen(rd());

    ImageRGB<T> im(static_cast<int>(im_size(0)),static_cast<int>(im_size(1)));
    im.parallel_for_all_pixels([&](ImageRGB<T>::PixelType &p,int i,int j)
    {
        // x = vec(0) between [-center(0),center(0)]
        // y = vec(1) between [-center(1),center(1)]
        Vec2 vec =  borns * Vec2(T(i) / T(im.width()), T(j) / T(im.height())) - center;

        T f(0), f2(0);

        T ax = vec(0) - s(0)/T(2);
        T bx = vec(0) + s(0)/T(2);
        T ay = vec(1) - s(1)/T(2);
        T by = vec(1) + s(1)/T(2);
        std::uniform_real_distribution<T> dis_x(ax, bx);
        std::uniform_real_distribution<T> dis_y(ay, by);


        int nb_sample = 100;
        for (int u =0; u < nb_sample; ++u) {
                T x = dis_x(gen);
                T y = dis_y(gen);
                T r = n.basic2D(x,y) ;
                f += r;
                f2 += r * r;
        }
        f /= nb_sample;
        f2 /= nb_sample;
        T sigma = std::sqrt(f2 - f*f);
        p = cm.map(f,sigma); // our filtered
//        p = ImageRGB<T>::itkPixel(cm.map(f)); // naive filtered
//        p = ImageRGB<T>::itkPixel(cm.map(n.basic2D(vec(0),vec(1)))); // unfiltered
    });

    return im;
}

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

    cm.filter(512,512,100,T(0.5));

    std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "filtrage de la color map timing: " << elapsed_seconds.count() << " s." << std::endl;

    ImageRGB<T> c0_ = cm.get_filtered();
    IO::save01_in_u8(c0_,TEMPO_PATH + "color_map_filtered.png");

//    ImageRGB<T> c0_;
//    IO::loadu8_in_01(c0_,TEMPO_PATH+ "color_map_filtered.png");
//    cm.set_filtered(c0_,T(0.5));

    Vec2 w_size(2,2);
    Vec2 im_size(512,512);
//    ImageGray<T> im_f = computeNoiseIMG(w_size,im_size,noise);

//    //compute f2
//    ImageGray<T> im_f2(im_f.width(),im_f.height());
//    im_f2.parallel_for_all_pixels([&im_f](ImageGray<T>::PixelType &p,int i,int j){
//        ImageGray<T>::PixelType v = im_f.pixelAbsolute(i,j);
//        p = v * v;
//    });

//    Mipmap<ImageGray<T>> m(im_f);
//    Mipmap<ImageGray<T>> m2(im_f2);

//    m.generate();
//    m2.generate();

    auto start_chrono2 = std::chrono::system_clock::now();

    ImageRGB<T> filtered = computeNoiseIMG(w_size,im_size,noise,cm);

    std::chrono::duration<double> elapsed_seconds2 = std::chrono::system_clock::now() - start_chrono2;
    std::cout << "synthe timing: " << elapsed_seconds2.count() << " s." << std::endl;

//    ImageGray<T> fmp;
//    m.fullMipmap(fmp);

//    ImageGray<T> fmp2;
//    m2.fullMipmap(fmp2);


//    IO::save01_in_u8(fmp,TEMPO_PATH + "noise.png");
//    IO::save01_in_u8(fmp2,TEMPO_PATH + "noise2.png");

    IO::save01_in_u8(filtered,TEMPO_PATH + "noise_filtered.png");

    return 0;
}
