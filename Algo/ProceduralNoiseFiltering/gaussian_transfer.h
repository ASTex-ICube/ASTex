#ifndef GAUSSIAN_TRANSFER_H
#define GAUSSIAN_TRANSFER_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <ASTex/image_rgb.h>
#include <ASTex/image_gray.h>
#include <ASTex/utils.h>

namespace ASTex {

template <typename IMG>
class Gaussian_transfer
{
public:
    using DataType = typename IMG::DataType;
    using PixelType = typename IMG::PixelType;

private:
    using T = DataType;
    struct PixelSortStruct
    {
        int x;
        int y;
        DataType value;
        bool operator < (const PixelSortStruct& other) const
        {
            return (value < other.value);
        }
    };

    static T Erf(T x)
    {
        // Save the sign of x
        int sign = 1;
        if (x < 0)
            sign = -1;
        x = std::abs(x);

        // A&S formula 7.1.26
        T t = T(1) / (T(1) + T(0.3275911) * x);
        T y = T(1) - (((((T(1.061405429) * t + -T(1.453152027)) * t) + T(1.421413741))
            * t + -T(0.284496736)) * t + T(0.254829592)) * t * std::exp(-x * x);

        return sign * y;
    }

    static T ErfInv(T x)
    {
        T w, p;
        w = -std::log((T(1) - x) * (T(1) + x));
        if (w < T(5))
        {
            w = w - T(2.5);
            p = T(2.81022636e-08);
            p = T(3.43273939e-07) + p * w;
            p = -T(3.5233877e-06) + p * w;
            p = -T(4.39150654e-06) + p * w;
            p = T(0.00021858087) + p * w;
            p = -T(0.00125372503) + p * w;
            p = -T(0.00417768164) + p * w;
            p = T(0.246640727) + p * w;
            p = T(1.50140941) + p * w;
        }
        else
        {
            w = std::sqrt(w) - T(3);
            p = -T(0.000200214257);
            p = T(0.000100950558) + p * w;
            p = T(0.00134934322) + p * w;
            p = -T(0.00367342844) + p * w;
            p = T(0.00573950773) + p * w;
            p = -T(0.0076224613) + p * w;
            p = T(0.00943887047) + p * w;
            p = T(1.00167406) + p * w;
            p = T(2.83297682) + p * w;
        }
        return p * x;
    }

    static T CDF(T x, T mu, T sigma)
    {
        T U = T(0.5) * (T(1) + Erf((x-mu)/(sigma* std::sqrt(T(2)))));
        return U;
    }

    static T invCDF(T U, T mu, T sigma)
    {
        T x = sigma*std::sqrt(T(2)) * ErfInv(T(2)*U-T(1)) + mu;
        return x;
    }

public:
    Gaussian_transfer();

    static void ComputeTinput(const IMG &input, IMG &T_input)
    {
        DataType * pix;
        for (unsigned int channel = 0; channel < IMG::NB_CHANNELS ; channel++) {
            // Sort pixels of example image
            std::vector<PixelSortStruct> sortedInputValues;
            sortedInputValues.resize(input.width() * input.height());

            input.for_all_pixels([&](const PixelType &p,int x,int y){
                sortedInputValues[y * input.width() + x].x = x;
                sortedInputValues[y * input.width() + x].y = y;

                PixelType tmp = p;
                pix = reinterpret_cast<DataType*>(&tmp);

                sortedInputValues[y * input.width() + x].value = pix[channel];
            });
            std::sort(sortedInputValues.begin(), sortedInputValues.end());

            // Assign Gaussian value to each pixel
            for (unsigned int i = 0; i < sortedInputValues.size() ; i++)
            {
                // Pixel coordinates
                int x = sortedInputValues[i].x;
                int y = sortedInputValues[i].y;
                // Input quantile (given by its order in the sorting)
                T U = (i + T(0.5)) / (sortedInputValues.size());
                // Gaussian quantile
                T G = invCDF(U, T(0.5), T(1)/T(6));
                G = clamp_scalar(G,T(0),T(1));
                // Store
                PixelType &p = T_input.pixelAbsolute(x, y);
                pix = reinterpret_cast<DataType*>(&p);
                pix[channel] = G;
            }
        }
    }

    static void ComputeinvT(const IMG& input, IMG & Tinv)
    {
        DataType * pix;
        for (unsigned int channel = 0; channel < IMG::NB_CHANNELS ; channel++) {
            // Sort pixels of example image
            std::vector<DataType> sortedInputValues;
            sortedInputValues.resize(input.width() * input.height());

            input.for_all_pixels([&](const PixelType &p, int x, int y)
            {
                PixelType tmp = p;
                pix = reinterpret_cast<DataType*>(&tmp);
                sortedInputValues[y * input.width() + x] = pix[channel];
            });
            std::sort(sortedInputValues.begin(), sortedInputValues.end());

            // Generate Tinv look-up table
            for (int i = 0; i < Tinv.width(); i++)
            {
                // Gaussian value in [0, 1]
                T G = (i + T(0.5)) / (Tinv.width());
                // Quantile value
                T U = CDF(G, T(0.5), T(1)/T(6));
                // Find quantile in sorted pixel values
                int index = int(std::floor(U * sortedInputValues.size()));
                // Get input value
                DataType value = sortedInputValues[index];
                // Store in LUT
                PixelType &p = Tinv.pixelAbsolute(i, 0);
                pix = reinterpret_cast<DataType*>(&p);
                pix[channel] = value;
            }
        }
    }

    static void invT(PixelType &p, const IMG &Tinv)
    {
        p = clamp_scalar(p,T(0),T(1));
        int size = Tinv.width() - 1;
        DataType *pix = reinterpret_cast<DataType*>(&p);
        DataType *pix_Tinv;
        for (unsigned int channel = 0; channel < IMG::NB_CHANNELS; channel++) {
            DataType value = pix[channel] * size;
//            DataType value_floor = std::floor(value);
//            DataType t = value - value_floor;

            int index = int(std::round(value));

            PixelType tmp = Tinv.pixelAbsolute(index, 0);
            pix_Tinv = reinterpret_cast<DataType*>(&tmp);

            pix[channel] = pix_Tinv[channel];

        }
    }

    static IMG invT(const IMG &input, const IMG &Tinv)
    {
        IMG output(input.width(), input.height());
        output.parallel_for_all_pixels([&](PixelType &pix,int x,int y){
            pix = input.pixelAbsolute(x,y);
            invT(pix,Tinv);
        });
        return output;
    }

};

}

#endif // GAUSSIAN_TRANSFER_H
