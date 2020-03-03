#include "gaussian_transfer.h"
#include <ASTex/utils.h>

namespace ASTex {


Gaussian_transfer::Gaussian_transfer()
{

}

float Gaussian_transfer::Erf(float x)
{
    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = std::abs(x);

    // A&S formula 7.1.26
    float t = 1.0f / (1.0f + 0.3275911f * x);
    float y = 1.0f - (((((1.061405429f * t + -1.453152027f) * t) + 1.421413741f)
        * t + -0.284496736f) * t + 0.254829592f) * t * std::exp(-x * x);

    return sign * y;
}

float Gaussian_transfer::ErfInv(float x)
{
    float w, p;
    w = -std::log((1.0f - x) * (1.0f + x));
    if (w < 5.000000f)
    {
        w = w - 2.500000f;
        p = 2.81022636e-08f;
        p = 3.43273939e-07f + p * w;
        p = -3.5233877e-06f + p * w;
        p = -4.39150654e-06f + p * w;
        p = 0.00021858087f + p * w;
        p = -0.00125372503f + p * w;
        p = -0.00417768164f + p * w;
        p = 0.246640727f + p * w;
        p = 1.50140941f + p * w;
    }
    else
    {
        w = std::sqrt(w) - 3.000000f;
        p = -0.000200214257f;
        p = 0.000100950558f + p * w;
        p = 0.00134934322f + p * w;
        p = -0.00367342844f + p * w;
        p = 0.00573950773f + p * w;
        p = -0.0076224613f + p * w;
        p = 0.00943887047f + p * w;
        p = 1.00167406f + p * w;
        p = 2.83297682f + p * w;
    }
    return p * x;
}

 float Gaussian_transfer::CDF(float x, float mu, float sigma)
{
    float U = 0.5f * (1 + Erf((x-mu)/(sigma* std::sqrt(2.0f))));
    return U;
}

float Gaussian_transfer::invCDF(float U, float mu, float sigma)
{
    float x = sigma*sqrtf(2.0f) * ErfInv(2.0f*U-1.0f) + mu;
    return x;
}

 void Gaussian_transfer::ComputeTinput(ImageRGBf &input, ImageRGBf &T_input)
{
    for (int channel = 0; channel < 3; channel++) {
        // Sort pixels of example image
        std::vector<PixelSortStruct> sortedInputValues;
        sortedInputValues.resize(input.width() * input.height());
        for (int y = 0; y < input.height(); y++)
        for (int x = 0; x < input.width(); x++)
        {
            sortedInputValues[y * input.width() + x].x = x;
            sortedInputValues[y * input.width() + x].y = y;
            sortedInputValues[y * input.width() + x].value = input.pixelAbsolute(x,y).GetNthComponent(channel);
        }
        std::sort(sortedInputValues.begin(), sortedInputValues.end());

        // Assign Gaussian value to each pixel
        for (unsigned int i = 0; i < sortedInputValues.size() ; i++)
        {
            // Pixel coordinates
            int x = sortedInputValues[i].x;
            int y = sortedInputValues[i].y;
            // Input quantile (given by its order in the sorting)
            float U = (i + 0.5f) / (sortedInputValues.size());
            // Gaussian quantile
            float G = invCDF(U, 0.5f, 0.16666f);
            G = clamp_scalar(G,0,1);
            // Store
            T_input.pixelAbsolute(x,y).SetNthComponent(channel,G);
        }
    }
}

void Gaussian_transfer::ComputeTinput(ImageGrayf &input, ImageGrayf &T_input){
        // Sort pixels of example image
        std::vector<PixelSortStruct> sortedInputValues;
        sortedInputValues.resize(input.width() * input.height());
        for (int y = 0; y < input.height(); y++)
        for (int x = 0; x < input.width(); x++)
        {
            sortedInputValues[y * input.width() + x].x = x;
            sortedInputValues[y * input.width() + x].y = y;
            sortedInputValues[y * input.width() + x].value = input.pixelAbsolute(x,y);
        }
        std::sort(sortedInputValues.begin(), sortedInputValues.end());

        // Assign Gaussian value to each pixel
        for (unsigned int i = 0; i < sortedInputValues.size() ; i++)
        {
            // Pixel coordinates
            int x = sortedInputValues[i].x;
            int y = sortedInputValues[i].y;
            //if( x == 126 && y==68)
            //{
            // Input quantile (given by its order in the sorting)
            float U = (i + 0.5f) / (sortedInputValues.size());
            // Gaussian quantile
            float G = invCDF(U, 0.5f, 0.16666f);

            G = clamp_scalar(G,0,1);

            // Store
            T_input.pixelAbsolute(x,y) = G;
            //}
        }
}

void Gaussian_transfer::ComputeinvT(ImageGrayf& input, ImageGrayf &Tinv)
{
    // Sort pixels of example image
    std::vector<float> sortedInputValues;
    sortedInputValues.resize(input.width() * input.height());
    for (int y = 0; y < input.height(); y++)
    for (int x = 0; x < input.width(); x++)
    {
        sortedInputValues[y * input.width() + x] = input.pixelAbsolute(x, y);
    }
    std::sort(sortedInputValues.begin(), sortedInputValues.end());

    // Generate Tinv look-up table
    for (int i = 0; i < Tinv.width(); i++)
    {
        // Gaussian value in [0, 1]
        float G = (i + 0.5f) / (Tinv.width());
        // Quantile value
        float U = CDF(G, 0.5f, 0.16666f);
        // Find quantile in sorted pixel values
        int index = (int)std::floor(U * sortedInputValues.size());
        // Get input value
        float I = sortedInputValues[index];
        // Store in LUT
        Tinv.pixelAbsolute(i, 0) = I;
    }
}

void Gaussian_transfer::ComputeinvT(ImageRGBf& input, ImageRGBf &Tinv)
{
    for (int channel = 0; channel < 3 ; channel++) {
        // Sort pixels of example image
        std::vector<float> sortedInputValues;
        sortedInputValues.resize(input.width() * input.height());
        for (int y = 0; y < input.height(); y++)
        for (int x = 0; x < input.width(); x++)
        {
            sortedInputValues[y * input.width() + x] = input.pixelAbsolute(x, y).GetNthComponent(channel);
        }
        std::sort(sortedInputValues.begin(), sortedInputValues.end());

        // Generate Tinv look-up table
        for (int i = 0; i < Tinv.width(); i++)
        {
            // Gaussian value in [0, 1]
            float G = (i + 0.5f) / (Tinv.width());
            // Quantile value
            float U = CDF(G, 0.5f, 0.16666f);
            // Find quantile in sorted pixel values
            int index = (int)std::floor(U * sortedInputValues.size());
            // Get input value
            float I = sortedInputValues[index];
            // Store in LUT
            Tinv.pixelAbsolute(i, 0).SetNthComponent(channel, I);
        }
    }
}

} //end namespace
