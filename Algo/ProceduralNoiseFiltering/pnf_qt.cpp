#include <ASTex/noise.h>
#include <ASTex/color_map.h>
#include <Eigen/Eigen>
#include <ASTex/mipmap.h>
#include "imageviewer.h"

using namespace ASTex;
using T = double;
using Vec2 = Eigen::Matrix<T,2,1>;
using Mat22 = Eigen::Matrix<T,2,2>;
using Color = Color_map<T>::Color;

inline ImageGray<T> computeNoiseIMG(const T &size_x,const T &size_y,const Noise<T> &n){
    Mat22 borns;
    borns << size_x, 0, 0, size_y;
    Vec2 center = Vec2(borns(0),borns(3)) / T(2);

    ImageGray<T> im(512,512);
    im.parallel_for_all_pixels([&](ImageGray<T>::PixelType &p,int i,int j)
    {
        // x = vec(0) between [-center(0),center(0)]
        // y = vec(1) between [-center(1),center(1)]
        Vec2 vec =  borns * Vec2(T(i) / T(im.width()), T(j) / T(im.height())) - center;

        // noise value at (x,y), vec(0) =x et vec(1)=y
        T r = n.basic2D(vec(0),vec(1));
        p = r;
        //im_f2.pixelAbsolute(i,j) = r * r;
    });

    return im;
}

int main(int argc, char **argv)
{
    QApplication app(argc, argv);

    Noise<T> noise;
    T size_x(6);
    T size_y(6);
    ImageGray<T> f = computeNoiseIMG(size_x,size_y,noise);

    auto iv = image_viewer(f,"noise",&app);

    iv->set_mouse_cb([&f](int b, int x, int y)
    {
        std::cout << "noise value = "<< f.pixelAbsolute(x,y) << " at " << x << " , " << y <<std::endl;
    });

    iv->set_key_cb([&](int code, char c)
    {
        std::ignore = code; // for warning
        switch(c)
        {
        case '+' :
            if(size_x > 0 && size_y > 0){
                size_x--;
                size_y--;
                f = computeNoiseIMG(size_x,size_y,noise);
                iv->update(f);
                iv->show();
            }
            break;
        case '-' :
                size_x++;
                size_y++;
                f = computeNoiseIMG(size_x,size_y,noise);
                iv->update(f);
                iv->show();
            break;
        default:
            std::cout << "key "<< c << std::endl;
        }
    });

    return app.exec();
}
