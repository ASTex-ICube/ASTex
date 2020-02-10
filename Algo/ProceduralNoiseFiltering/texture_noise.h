#ifndef TEXTURE_NOISE_H
#define TEXTURE_NOISE_H

#include <string>
#include <ASTex/easy_io.h>
#include <ASTex/image_spectral.h>
#include <ASTex/rpn_utils.h>
#include <ASTex/color_map.h>
#include <ASTex/exr_io.h>

namespace ASTex {


template <typename T>
class TextureNoise
{
private:
    ImageGray<T> noise;
public:
    using Vec2 = Eigen::Matrix<T,2,1>;
    using Color = typename Color_map<T>::Color;

    TextureNoise() {}
    TextureNoise(ImageSpectrald &psd) {
        noise.initItk(psd.width(),psd.height());

        ImageSpectrald phase;
        rpn_scalar(psd, phase, noise);

        noise.for_all_pixels([] (typename ImageGray<T>::PixelType &pix)
        {
            pix = clamp_scalar(pix, T(0), T(1));
        });
    }

    TextureNoise(ImageSpectrald &psd, const T &m, const T &s) {

//        T var(0),mean(0);
//        psd.for_all_pixels([&](typename ImageGray<T>::PixelType &pix, int x, int y){
//            if( x == psd.width() /2 && y == psd.height()/2)
//                mean = pix;
//            else
//                var += pix * pix ;
//        });

        noise.initItk(psd.width(),psd.height());

        ImageSpectrald phase;
        rpn_scalar(psd, phase, noise);

        //IO::save_phase(phase, TEMPO_PATH + "phases.png");
        IO::save01_in_u8(psd, TEMPO_PATH + "psd.png");
        //IO::EXR::save(psd, TEMPO_PATH + "psd.exr");
        T mean = getMean(noise);
        T sigma = getStDev(noise);
        std::cout << mean << std::endl;
        std::cout << sigma << std::endl;

        noise.for_all_pixels([&] (typename ImageGray<T>::PixelType &pix)
        {
            pix *= s * 1.0/sigma;
            pix += m - mean;
            pix = clamp_scalar(pix, T(0), T(1));
        });

        std::cout << getMean(noise) << std::endl;
        std::cout << getStDev(noise) << std::endl;

    }

    ImageGray<T> getNoise() const {
        return noise;
    }

    void setNoise(const ImageGray<T> &im){
        noise = im;
    }

    int width() const {
        return noise.width();
    }

    int height() const {
        return noise.height();
    }

    T get(const int &i, const int &j) const {
        int x = (i % width() + width()) % width();
        int y = (j % height() + height()) % height();

        return noise.pixelAbsolute(x, y);
    }

    T get(const Vec2 &pos) const {
        int x = static_cast<int>(pos(0));
        int y = static_cast<int>(pos(1));
        return get(x,y);
    }

private :
    template<typename T2,typename init, typename func>
    T2 get_over_footprint(const Vec2 &pos, const Vec2 &footprint, const init &initialisation, const func &f) const{
        Vec2 corner = pos - footprint * T(0.5);
        T2 ret = initialisation();
        T nb(0);
        for(int i = 0; i <= static_cast<int>(footprint(1)) ; ++i){
            int y = static_cast<int>(corner(1)) + i;
            for (int j=0; j <= static_cast<int>(footprint(0)); ++j) {
                int x = static_cast<int>(corner(0)) + j;
                T value_noise = get(x,y);
                ret += f(value_noise);
                nb++;
            }
        }

        ret /= nb;

        return ret;
    }

public:
    T get_noise_mean_over_footprint(const Vec2 &pos, const Vec2 &footprint) const{
        return get_over_footprint<T>(pos, footprint, [](){return T(0);}, [](const T&x) {
            return x;
        });
    }

    T get_squared_noise_mean_over_footprint(const Vec2 &pos, const Vec2 &footprint) const{
        return get_over_footprint<T>(pos, footprint, [](){return T(0);}, [](const T&x) {
            return x * x;
        });
    }

    Color get_color_mean_over_footprint(const Vec2 &pos, const Vec2 &footprint, const Color_map<T> &cm) const{
        return get_over_footprint<Color>(pos, footprint, [](){return Color(0,0,0);}, [&](const T&x) {
            return cm.map(x, 0);
        });
    }
};

}
#endif // TEXTURE_NOISE_H
