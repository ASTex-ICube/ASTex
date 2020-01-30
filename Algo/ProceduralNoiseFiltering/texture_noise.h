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
    TextureNoise() {}
    TextureNoise(const std::string &filename) {
        ImageSpectrald psd;
        IO::loadu8_in_01(psd, filename);

        noise.initItk(psd.width(),psd.height());
        psd.for_all_pixels([] (ImageSpectrald::PixelType &pix)
        {
            if(pix>0)
                pix /= pix;
        });

        ImageSpectrald phase;
        rpn_scalar(psd, phase, noise);
        noise.for_all_pixels([] (typename ImageGray<T>::PixelType &pix)
        {
            pix = clamp_scalar(pix, T(0), T(1));
        });
    }
};

}
#endif // TEXTURE_NOISE_H
