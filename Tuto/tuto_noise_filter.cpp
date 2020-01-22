#include <ASTex/noise.h>
#include <ASTex/easy_io.h>
#include <Eigen/Eigen>

using namespace ASTex;

int main()
{
    using T = double;
    using Vec2 = Eigen::Matrix<T,2,1>;
    using Mat22 = Eigen::Matrix<T,2,2>;
    ImageGray<T> im(512,512);
    Noise<T> noise;

    im.parallel_for_all_pixels([&](T &p,int i,int j)
    {
        Mat22 borns;
        borns << 2, 0, 0, 2;
        Vec2 center = Vec2(borns(0),borns(3)) / T(2);
        Vec2 vec =  borns * Vec2(T(i) / T(im.width()), T(j) / T(im.height())) - center;
//        T x = T(i) / T(im.width()) * borns(0) - center(0);      // x between [-center,center]
//        T y = T(j) / T(im.height()) * borns(3) - center(1);     // y between [-center,center]
        T r = noise.basic2D(vec(0),vec(1));                               // noise value at (x,y)
        p = r;
    });

    IO::save01_in_u8(im,TEMPO_PATH + "noise.png");

    return 0;
}
