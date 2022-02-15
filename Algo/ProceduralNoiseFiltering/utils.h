#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <ASTex/easy_io.h>

namespace ASTex {

template <typename IMG>
void gammaCorr(IMG &img)
{
    typename IMG::DataType *p;
    img.for_all_pixels([&](typename IMG::PixelType & pix){
        p = reinterpret_cast<typename IMG::DataType*>(&pix);
        for (unsigned int channel = 0; channel < IMG::NB_CHANNELS ; channel++)
            p[channel] = std::pow(p[channel], 0.45);
    });
}

template <typename IMG>
void gammaInvCorr(IMG &img)
{
    typename IMG::DataType *p;
    img.for_all_pixels([&](typename IMG::PixelType & pix){
        p = reinterpret_cast<typename IMG::DataType*>(&pix);
        for (unsigned int channel = 0; channel < IMG::NB_CHANNELS ; channel++)
            p[channel] = std::pow(p[channel], 2.2);
    });
}

template <typename IMG>
void loadGamma(IMG & img, const std::string &filename ,const bool & _16 = false)
{
    if(_16)
        IO::loadu16_in_01(img, filename);
    else
        IO::loadu8_in_01(img, filename);

    gammaInvCorr(img);
}

template <typename IMG>
void saveGamma(IMG & img, const std::string &filename ,const bool & _16 = false)
{
    gammaCorr(img);

    if(_16)
        IO::save01_in_u16(img, filename);
    else
        IO::save01_in_u8(img, filename);
}

}

#endif // UTILS_H
