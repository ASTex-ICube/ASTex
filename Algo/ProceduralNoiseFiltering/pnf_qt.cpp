#include <ASTex/noise.h>
#include <ASTex/color_map.h>
#include <Eigen/Eigen>
#include <ASTex/rpn_utils.h>
#include "imageviewer.h"

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
    std::mt19937 gen(0);

    ImageRGB<T> im(static_cast<int>(im_size(0)),static_cast<int>(im_size(1)));
    im.for_all_pixels([&](ImageRGB<T>::PixelType &p,int i,int j)
    {
        // x = vec(0) between [-center(0),center(0)]
        // y = vec(1) between [-center(1),center(1)]
        Vec2 vec =  borns * Vec2(T(i) / T(im.width()), T(j) / T(im.height())) - center;

        T f(0), f2(0);

        if(i==41 && j==346)
            std::cerr << "toto" << std::endl;

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
        p = cm.map(f,sigma);
    });

    return im;
}

int main(int argc, char **argv)
{
    QApplication app(argc, argv);

//    Noise<T> noise;
//    T size_x(6);
//    T size_y(6);
//    ImageGray<T> psd;
//    psd.load(TEMPO_PATH + "spectra/donut_black.png");

//    ImageGray<T> f(512,512) ;
//    gray_RPN(psd,f,0,0);
//    IO::save01_in_u8(f,TEMPO_PATH + "noise_psd.png");

    Noise<T> noise;
    Color_map<T> cm;
    cm.add_color(0,Color(1,1,0));
    cm.add_color(1,Color(1,0,0));
    cm.add_color(2,Color(0,0,0));
    cm.add_color(3,Color(1,1,1));

    ImageRGB<T> c0_;
    IO::loadu8_in_01(c0_,TEMPO_PATH+ "color_map_filtered.png");
    cm.set_filtered(c0_,T(0.5));

    Vec2 w_size(2,2);
    Vec2 im_size(512,512);

    ImageRGB<T> f = computeNoiseIMG(w_size,im_size,noise,cm);

    auto iv = image_viewer(f,"noise",&app);

    iv->set_mouse_cb([&f](int b, int x, int y)
    {
        std::cout << "color value = "<< f.pixelAbsolute(x,y) << " at " << x << " , " << y <<std::endl;
    });

    iv->set_key_cb([&](int code, char c)
    {
        std::ignore = code; // for warning
        switch(c)
        {
        case '+' :
            if(w_size(0) > 0 && w_size(1) > 0){
                w_size(0)--;
                w_size(1)--;
                f = computeNoiseIMG(w_size,im_size,noise,cm);
                iv->update(f);
                iv->show();
            }
            break;
        case '-' :
                w_size(0)++;
                w_size(1)++;
                f = computeNoiseIMG(w_size,im_size,noise,cm);
                iv->update(f);
                iv->show();
            break;
        default:
            std::cout << "key "<< c << std::endl;
        }
    });

    return app.exec();
}
