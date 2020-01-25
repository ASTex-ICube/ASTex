#include <ASTex/noise.h>
#include <ASTex/easy_io.h>
#include <ASTex/color_map.h>
#include <Eigen/Eigen>
#include <ASTex/mipmap.h>

using namespace ASTex;

int main(int argc, char **argv)
{
    //QApplication app(argc, argv);
    using T = double;
    using Vec2 = Eigen::Matrix<T,2,1>;
    using Mat22 = Eigen::Matrix<T,2,2>;
    using Color = Color_map<T>::Color;
    ImageGray<T> im_f(512,512);
    ImageGray<T> im_f2(512,512);

    ImageRGB<T> filtered(512,512);

    Noise<T> noise;
    Color_map<T> cm;
    cm.add_color(0,Color(0,0,0));
    cm.add_color(1,Color(0,0,1));
    cm.add_color(3,Color(0,1,0));
    cm.add_color(4,Color(1,0,0));
    cm.add_color(6,Color(1,1,1));

    //cm.export_palette(TEMPO_PATH + "palette.gnu");
    //cm.export_courbe(TEMPO_PATH + "data.txt");

    ImageRGB<T> c0_ = cm.filter(512,512,10,T(1)/T(2));

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

    //auto iv = image_view(f,"noise",&app);

    return 0;
    //return app.exec();
}
