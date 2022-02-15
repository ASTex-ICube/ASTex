#ifndef _AUTOREGRESSIVE_H_
#define _AUTOREGRESSIVE_H_

#include "histogram.h"
#include "easy_io.h"

template<typename I>
class Autoregressive
{
public:
	using ImageType					= I;
	using HitMapType				= ImageGrayu8;
	using WhiteNoiseImageType		= I;
	using WhiteNoisePixelType		= typename I::PixelType;
	using PixelType					= typename ImageType::PixelType;
	using DataType					= typename ImageType::DataType;
	using PixelPosType				= itk::Index<2>;

	Autoregressive();

	ImageType simulateFromLeftUp() const;
	void setOrder(int p, int q = 0);
	void setCoeff(int i, int j, const PixelType &value);
	void setCoeff(std::function<PixelType (unsigned int, unsigned int)> functor);
	PixelType coeff(int i, int j) const;

	void setWhiteNoiseParameters(WhiteNoisePixelType mean, WhiteNoisePixelType variance);
	void setWhiteNoiseMean(WhiteNoisePixelType mean);
	void setWhiteNoiseVariance(WhiteNoisePixelType variance);
	void setConstant(const PixelType &value);
	void setWidth(unsigned int width);
	void setHeight(unsigned int height);

private:

	void compute(ImageType &image, int x, int y, const WhiteNoiseImageType &whiteNoise, const HitMapType *hitMap) const;

	ImageType m_coeffs;
	int m_width;
	int m_height;
	WhiteNoisePixelType m_wnMean;
	WhiteNoisePixelType m_wnVariance;
	PixelType m_constant;
};

template<typename I>
Autoregressive<I>::Autoregressive() :
	m_coeffs(),
	m_width(0),
	m_height(0),
	m_wnMean(),
	m_wnVariance()
{}

template<typename I>
void Autoregressive<I>::compute(ImageType &image, int x, int y, const WhiteNoiseImageType &whiteNoise, const HitMapType *hitMap) const
{
	image.pixelAbsolute(x, y) = ImageType::zero();
	int pMax = (m_coeffs.width() - 1)/2;
	int qMax = (m_coeffs.height() - 1)/2;
	for(int p = -pMax; p<pMax; ++p)
	{
		for(int q = -qMax; q<qMax; ++q)
		{
			if( (p != 0 || q != 0)
			&&	x+p >= 0 && x+p < image.width()
			&&	y+q >= 0 && y+q < image.height()
			&&	(hitMap == nullptr || hitMap->pixelAbsolute(x+p, y+q) == 1))
			{
				image.pixelAbsolute(x, y) += image.pixelAbsolute(x+p, y+q) * coeff(p, q) + whiteNoise.pixelAbsolute(x, y);
			}
		}
	}
}

template<typename I>
typename Autoregressive<I>::ImageType Autoregressive<I>::simulateFromLeftUp() const
{
	assert(m_width != 0 && m_height != 0);
	ImageType image;
	WhiteNoiseImageType W = whiteNoise<WhiteNoiseImageType>(m_width, m_height, m_wnMean, m_wnVariance);
	HitMapType hitMap;
	hitMap.initItk(m_width, m_width, true);
	image.initItk(m_width, m_height, true);
	image.pixelAbsolute(0, 0) = m_constant;
	hitMap.pixelAbsolute(0, 0) = 1;
	int d = 1;
	for(d=1; d<m_width; ++d)
	{
		for(int x=0; x<=d; ++x)
		{
			compute(image, x, d, W, &hitMap);
			hitMap.pixelAbsolute(x, d) = 1;
		}
		for(int y=0; y<=d; ++y)
		{
			compute(image, d, y, W, &hitMap);
			hitMap.pixelAbsolute(d, y) = 1;
		}
	}
	for(int x=d+1; x<m_width; ++x)
	{
		for(int y=0; y<m_height; ++y)
		{
			compute(image, d, y, W, &hitMap);
			hitMap.pixelAbsolute(d, y) = 1;
		}
	}
	return image;
}

template<typename I>
void Autoregressive<I>::setOrder(int p, int q)
{
	assert(p > 0 && q >=0);
	if(q == 0)
		q = p;
	m_coeffs.initItk(p*2+1, q*2+1, true);
}

template<typename I>
void Autoregressive<I>::setCoeff(int p, int q, const PixelType &value)
{
	assert(p!=0 || q!= 0);
	unsigned int pMax = (m_coeffs.width() - 1)/2;
	unsigned int qMax = (m_coeffs.height() - 1)/2;
	m_coeffs.pixelAbsolute(p+pMax, q+qMax) = value;
}

template<typename I>
void Autoregressive<I>::setCoeff(std::function<PixelType (unsigned int, unsigned int)> functor)
{
	unsigned int pMax = (m_coeffs.width() - 1)/2;
	unsigned int qMax = (m_coeffs.height() - 1)/2;
	for(unsigned int i=0; i<m_coeffs.width(); ++i)
	{
		for(unsigned int j=0; j<m_coeffs.height(); ++j)
		{
			int p = i-pMax;
			int q = j-qMax;
			if(p != 0 || q != 0)
			{
				m_coeffs.pixelAbsolute(i, j) = ImageType::zero();
			}
			else
				m_coeffs.pixelAbsolute(i, j) = functor(p, q);
		}
	}
}

template<typename I>
typename Autoregressive<I>::PixelType Autoregressive<I>::coeff(int p, int q) const
{
	unsigned int pMax = (m_coeffs.width() - 1)/2;
	unsigned int qMax = (m_coeffs.height() - 1)/2;
	PixelType pix;
	return (pix = m_coeffs.pixelAbsolute(p+pMax, q+qMax));
}

template<typename I>
void Autoregressive<I>::setWhiteNoiseParameters(WhiteNoisePixelType mean, WhiteNoisePixelType variance)
{
	m_wnMean = mean;
	m_wnVariance = variance;
}

template<typename I>
void Autoregressive<I>::setWhiteNoiseMean(WhiteNoisePixelType mean)
{
	m_wnMean = mean;
}

template<typename I>
void Autoregressive<I>::setWhiteNoiseVariance(WhiteNoisePixelType variance)
{
	m_wnVariance = variance;
}

template<typename I>
void Autoregressive<I>::setConstant(const PixelType &value)
{
	m_constant = value;
}

template<typename I>
void Autoregressive<I>::setWidth(unsigned int width)
{
	m_width = width;
}

template<typename I>
void Autoregressive<I>::setHeight(unsigned int height)
{
	m_height = height;
}

#endif //_AUTOREGRESSIVE_H_
