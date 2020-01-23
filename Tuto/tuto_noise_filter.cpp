#include <ASTex/noise.h>
#include <ASTex/easy_io.h>
#include <ASTex/color_map.h>
#include <Eigen/Eigen>

using namespace ASTex;

int main()
{
    using T = double;
    using Vec2 = Eigen::Matrix<T,2,1>;
    using Mat22 = Eigen::Matrix<T,2,2>;
    ImageRGB<T> im(512,512);
    Noise<T> noise;
    Color_map<T> cm;
    cm.add_color(0,Color_map<T>::Color(0,0,0));
    cm.add_color(1,Color_map<T>::Color(0,0,1));
    cm.add_color(3,Color_map<T>::Color(0,1,0));
    cm.add_color(4,Color_map<T>::Color(1,0,0));
    cm.add_color(6,Color_map<T>::Color(1,1,1));

    im.for_all_pixels([&](ImageRGB<T>::PixelType &p,int i,int j)
    {
        Mat22 borns;
        borns << 20, 0, 0, 20;
        Vec2 center = Vec2(borns(0),borns(3)) / T(2);
        Vec2 vec =  borns * Vec2(T(i) / T(im.width()), T(j) / T(im.height())) - center;
//        T x = T(i) / T(im.width()) * borns(0) - center(0);      // x between [-center,center]
//        T y = T(j) / T(im.height()) * borns(3) - center(1);     // y between [-center,center]
        T r = noise.basic2D(vec(0),vec(1));                               // noise value at (x,y), vec(0) =x et vec(1)=y
        p = cm.map(T(i) / T(im.width()));
    });

    IO::save01_in_u8(im,TEMPO_PATH + "noise.png");

    return 0;
}
