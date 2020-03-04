#ifndef GAUSSIAN_TRANSFER_H
#define GAUSSIAN_TRANSFER_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <ASTex/image_rgb.h>
#include <ASTex/image_gray.h>
#include <ASTex/utils.h>
#include "ASTex/pca.h"

namespace ASTex {

template<typename I>
class Gaussian_transfer
{
public:

	using ImageType		= I;
	using PixelType		= typename ImageType::PixelType;
	using DataType		= typename ImageType::DataType;
	using LutType		= ImageType;

	Gaussian_transfer();
	static void ComputeTinput(ImageType &input, LutType &T_input);
	static void ComputeinvT(ImageType& input, LutType& Tinv);

private:
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

    static float Erf(float x);
    static float ErfInv(float x);

    static float CDF(float x, float mu, float sigma);
    static float invCDF(float U, float mu, float sigma);
};

template<typename I>
Gaussian_transfer<I>::Gaussian_transfer()
{}

template<typename I>
float Gaussian_transfer<I>::Erf(float x)
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

template<typename I>
float Gaussian_transfer<I>::ErfInv(float x)
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

template<typename I>
float Gaussian_transfer<I>::CDF(float x, float mu, float sigma)
{
	float U = 0.5f * (1 + Erf((x-mu)/(sigma* std::sqrt(2.0f))));
	return U;
}

template<typename I>
float Gaussian_transfer<I>::invCDF(float U, float mu, float sigma)
{
	float x = sigma*sqrtf(2.0f) * ErfInv(2.0f*U-1.0f) + mu;
	return x;
}

template<typename I>
void Gaussian_transfer<I>::ComputeTinput(ImageType &input, LutType &T_input)
{
	unsigned pixelSize = sizeof(PixelType)/sizeof(DataType);
	DataType *pixData;
	for (unsigned channel = 0; channel < pixelSize; channel++) {
		// Sort pixels of example image
		std::vector<PixelSortStruct> sortedInputValues;
		sortedInputValues.resize(input.width() * input.height());
		input.for_all_pixels([&] (PixelType &pix, int x, int y)
		{
			sortedInputValues[y * input.width() + x].x = x;
			sortedInputValues[y * input.width() + x].y = y;
			pixData = reinterpret_cast<DataType *>(&pix);
			sortedInputValues[y * input.width() + x].value = pixData[channel];
		});
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
			PixelType &pix = T_input.pixelAbsolute(x, y);
			pixData = reinterpret_cast<DataType *>(&pix);
			pixData[channel] = G;
		}
	}
}

template<typename I>
void Gaussian_transfer<I>::ComputeinvT(ImageType& input, LutType &Tinv)
{
	unsigned pixelSize = sizeof(PixelType)/sizeof(DataType);
	DataType *pixData;
	for (unsigned channel = 0; channel < pixelSize; channel++) {
		// Sort pixels of example image
		std::vector<DataType> sortedInputValues;
		sortedInputValues.resize(input.width() * input.height());
		input.for_all_pixels([&] (PixelType &pix, int x, int y)
		{
			pixData = reinterpret_cast<DataType *>(&pix);
			sortedInputValues[y * input.width() + x]= pixData[channel];
		});
		std::sort(sortedInputValues.begin(), sortedInputValues.end());

		// Generate Tinv look-up table
		for (int i = 0; i < Tinv.width(); i++)
		{
			// Gaussian value in [0, 1]
			float G = (i + 0.5f) / (Tinv.width());
			// Quantile value
			float U = CDF(G, 0.5f, 0.16666f);
			// Find quantile in sorted pixel values
			int index = int(std::floor(U * sortedInputValues.size()));
			// Get input value
			DataType value = sortedInputValues[index];
			// Store in LUT
			PixelType &pix = Tinv.pixelAbsolute(i, 0);
			pixData = reinterpret_cast<DataType *>(&pix);
			pixData[channel] = value;
		}
	}
}


}

#endif // GAUSSIAN_TRANSFER_H
