#ifndef GAUSSIAN_TRANSFER_H
#define GAUSSIAN_TRANSFER_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <ASTex/image_rgb.h>
#include <ASTex/image_gray.h>

namespace ASTex {

class Gaussian_transfer
{
private:
    struct PixelSortStruct
    {
        int x;
        int y;
        float value;
        bool operator < (const PixelSortStruct& other) const
        {
            return (value < other.value);
        }
    };

    static float Erf(float x);
    static float ErfInv(float x);

    static float CDF(float x, float mu, float sigma);
    static float invCDF(float U, float mu, float sigma);

public:
    Gaussian_transfer();
    static void ComputeTinput(ImageGrayf &input, ImageGrayf &T_input);
    static void ComputeTinput(ImageRGBf &input, ImageRGBf &T_input);

    static void ComputeinvT(ImageGrayf& input, ImageGrayf& Tinv);
    static void ComputeinvT(ImageRGBf& input, ImageRGBf& Tinv);
};

}

#endif // GAUSSIAN_TRANSFER_H
