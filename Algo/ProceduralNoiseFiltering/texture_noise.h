#ifndef TEXTURE_NOISE_H
#define TEXTURE_NOISE_H

#include <string>
#include <ASTex/easy_io.h>
#include <ASTex/image_spectral.h>
#include <ASTex/rpn_utils.h>

namespace ASTex {


template <typename T>
class TextureNoise
{
private:
    ImageGray<T> noise;
public:
    using Vec2 = Eigen::Matrix<T,2,1>;

    TextureNoise() {}
    TextureNoise(const ImageSpectrald psd) {

        ImageSpectrald sqrt_psd(psd.width(),psd.height());
        sqrt_psd.for_all_pixels([&] (ImageSpectrald::PixelType &pix,int i, int j)
        {
            ImageSpectrald::PixelType p = psd.pixelAbsolute(i,j);
            pix = std::sqrt(p);
        });

        ImageSpectrald phase;
        noise.initItk(psd.width(),psd.height());
        rpn_scalar(sqrt_psd, phase, noise);
        noise.for_all_pixels([] (typename ImageGray<T>::PixelType &pix)
        {
            pix = clamp_scalar(pix, T(0), T(1));
        });
        std::cout << getMean(noise) << std::endl;
        std::cout << getStDev(noise) << std::endl;
    }

    ImageGray<T> getNoise() const {
        return noise;
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
};

}
#endif // TEXTURE_NOISE_H
