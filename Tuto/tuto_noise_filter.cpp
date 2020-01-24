#include <ASTex/noise.h>
#include <ASTex/easy_io.h>
#include <ASTex/color_map.h>
#include <Eigen/Eigen>
#include <ASTex/mipmap.h>

using namespace ASTex;

int main()
{
    using T = double;
    using Vec2 = Eigen::Matrix<T,2,1>;
    using Mat22 = Eigen::Matrix<T,2,2>;
    ImageGray<T> im_f(512,512);
    ImageGray<T> im_f2(512,512);

    ImageRGB<T> filtered(512,512);

    Noise<T> noise;
    Color_map<T> cm;
    cm.add_color(0,Color_map<T>::Color(0,0,0));
    cm.add_color(1,Color_map<T>::Color(0,0,1));
    cm.add_color(3,Color_map<T>::Color(0,1,0));
    cm.add_color(4,Color_map<T>::Color(1,0,0));
    cm.add_color(6,Color_map<T>::Color(1,1,1));

    cm.export_palette(TEMPO_PATH + "palette.gnu");
    cm.export_courbe(TEMPO_PATH + "data.txt");


    ImageRGB<T> c0_(512,512);
    T sigma_max(6);
    c0_.parallel_for_all_pixels([&](ImageRGB<T>::PixelType &p,int i,int j)
    {
        T f = j / c0_.height();
        T sigma = i / c0_.width() * sigma_max;


    });

    Mat22 borns;
    borns << 2, 0, 0, 2;
    Vec2 center = Vec2(borns(0),borns(3)) / T(2);

    im_f.parallel_for_all_pixels([&](ImageGray<T>::PixelType &p,int i,int j)
    {
        // x = vec(0) between [-center(0),center(0)]
        // y = vec(1) between [-center(1),center(1)]
        Vec2 vec =  borns * Vec2(T(i) / T(im_f.width()), T(j) / T(im_f.height())) - center;

        // noise value at (x,y), vec(0) =x et vec(1)=y
        T r = noise.basic2D(vec(0),vec(1));
        p = r;
        im_f2.pixelAbsolute(i,j) = r * r;
    });

    Mipmap<ImageGray<T>> m(im_f);
    Mipmap<ImageGray<T>> m2(im_f2);

    m.generate();
    m2.generate();

    filtered.parallel_for_all_pixels([&] (ImageRGB<T>::PixelType &p, int i, int j)
    {
        p = cm.map(im_f.pixelAbsolute(i,j));
    });

    ImageGray<T> fmp;
    m.fullMipmap(fmp);

    ImageGray<T> fmp2;
    m2.fullMipmap(fmp2);


    IO::save01_in_u8(fmp,TEMPO_PATH + "noise.png");
    IO::save01_in_u8(fmp2,TEMPO_PATH + "noise2.png");
    IO::save01_in_u8(filtered,TEMPO_PATH + "noise_filtered.png");

    return 0;
}
