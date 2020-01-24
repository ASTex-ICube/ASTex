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

    cm.export_palette(TEMPO_PATH + "palette.gnu");

    im.parallel_for_all_pixels([&](ImageRGB<T>::PixelType &p,int i,int j)
    {
        Mat22 borns;
        borns << 2, 0, 0, 2;
        Vec2 center = Vec2(borns(0),borns(3)) / T(2);

        // x = vec(0) between [-center(0),center(0)]
        // y = vec(1) between [-center(1),center(1)]
        Vec2 vec =  borns * Vec2(T(i) / T(im.width()), T(j) / T(im.height())) - center;

        // noise value at (x,y), vec(0) =x et vec(1)=y
        T r = noise.basic2D(vec(0),vec(1));
        p = cm.map(r);
    });

    IO::save01_in_u8(im,TEMPO_PATH + "noise.png");

    return 0;
}
