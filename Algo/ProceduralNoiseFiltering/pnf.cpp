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

int main()
{
    Noise<T> noise;
    Color_map<T> cm;
    cm.add_color(0,Color(1,1,0));
    cm.add_color(1,Color(1,0,0));
    cm.add_color(2,Color(0,0,0));
    cm.add_color(3,Color(1,1,1));

    //cm.export_palette(TEMPO_PATH + "palette.gnu");
    cm.export_courbe(TEMPO_PATH + "data.txt");

    cm.filter(512,512,10,T(1)/T(2));

    ImageRGB<T> c0_ = cm.get_filtered();

    ImageGray<T> im_f = computeNoiseIMG(2,2,noise);

    ImageGray<T> im_f2(im_f.width(),im_f.height());
    im_f2.parallel_for_all_pixels([&im_f](ImageGray<T>::PixelType &p,int i,int j){
        ImageGray<T>::PixelType v = im_f.pixelAbsolute(i,j);
        p = v * v;
    });

    ImageRGB<T> filtered(im_f.width(),im_f.height());

    Mipmap<ImageGray<T>> m(im_f);
    Mipmap<ImageGray<T>> m2(im_f2);

    m.generate();
    m2.generate();

    auto start_chrono = std::chrono::system_clock::now();
    filtered.parallel_for_all_pixels([&] (ImageRGB<T>::PixelType &p, int i, int j)
    {
        p = ImageRGB<T>::itkPixel(cm.map(im_f.pixelAbsolute(i,j)));
    });
    std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
    std::cout << "synthe timing: " << elapsed_seconds.count() << " s." << std::endl;

    ImageGray<T> fmp;
    m.fullMipmap(fmp);

    ImageGray<T> fmp2;
    m2.fullMipmap(fmp2);


    IO::save01_in_u8(fmp,TEMPO_PATH + "noise.png");
    IO::save01_in_u8(fmp2,TEMPO_PATH + "noise2.png");
    IO::save01_in_u8(c0_,TEMPO_PATH + "color_map_filtered.png");
    IO::save01_in_u8(filtered,TEMPO_PATH + "noise_filtered.png");

    return 0;
}
