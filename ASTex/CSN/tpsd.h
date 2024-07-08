#ifndef TPSD_H_
#define TPSD_H_

#include <Eigen/Core>
#include "Algo/ProceduralNoiseFiltering/gaussian_transfer.h"
#include "ASTex/utils.h"
#include "ASTex/pca.h"
#include <cmath>
#include "ASTex/colorspace_filters.h"
#include "ASTex/easy_io.h"
#include "ASTex/histogram.h"
#include "ASTex/texture_pool.h"

//#define CSN_ENABLE_DEBUG
#ifdef CSN_ENABLE_DEBUG
#define print_debug(x) std::cout << x << std::endl
#define save_debug(x, y) IO::save01_in_u8(x, y)
#else
#define print_debug(x)
#define save_debug(x, y)
#endif

namespace ASTex
{

namespace CSN
{

template<typename I>
class TPSD
{
public:

	using ImageType					= I;
	using PixelType					= typename ImageType::PixelType;
	using DataType					= typename ImageType::DataType;
	using PixelPosType				= itk::Index<2>;
	using TPSDFunctionType			= std::function<std::function<void (PixelType&, int, int)> (double, double)>;

	//step 1
	void setNbSamples(unsigned nbSamplesX, unsigned nbSamplesY);
	void setSize(unsigned width, unsigned height);
	void setParallelogramCheat(bool useParallelogramCheat, bool isOnX=true);

	//step 2
	/**
	 * @brief setOperatorFunction provides a function describing how the T-PSD should be computed according to the periodic coordinates x and y.
	 * @param f is a function taking for argument a for_all_pixels type function; an x coordinate between 0 and 1; and a y coordinate between 0 and 1.
	 */
	void setOperatorFunction(TPSDFunctionType &f);

	/**
	 * @brief computeExplicitTPSD uses the provided function and the number of samples to generate an explicit TPSD.
	 * WARNING: DO NOT USE WITH TOO MANY SAMPLES: MAY CAUSE STACK SMASHES. KNOW YOUR AVAILABLE MEMORY AND THE SIZE OF YOUR SPECTRA.
	 */
	void computeExplicitTPSD();

	//step 3
	const ImageType operator()(unsigned x, unsigned y) const;

private:

	TPSDFunctionType m_function;
	TexturePool<ImageType, double> m_explicitTPSD;
	unsigned m_nbSamplesX, m_nbSamplesY;
	unsigned m_width, m_height;
	bool m_parallelogramCheat;
	bool m_parallelogramOnX;
};

template<typename I>
void TPSD<I>::setNbSamples(unsigned nbSamplesX, unsigned nbSamplesY)
{
	assert(nbSamplesX>1 && nbSamplesY>1);
	m_nbSamplesX = nbSamplesX;
	m_nbSamplesY = nbSamplesY;
}

template<typename I>
void TPSD<I>::setSize(unsigned width, unsigned height)
{
	assert(width>1 && height>1);
	m_width = width;
	m_height = height;
}

template<typename I>
void TPSD<I>::setParallelogramCheat(bool useParallelogramCheat, bool isOnX)
{
	m_parallelogramCheat = useParallelogramCheat;
	m_parallelogramOnX = isOnX;
}

template<typename I>
void TPSD<I>::setOperatorFunction(TPSDFunctionType &f)
{
	m_function = f;
}

template<typename I>
void TPSD<I>::computeExplicitTPSD()
{
	for(unsigned x=0; x<m_nbSamplesX; ++x)
		for(unsigned y=0; y<m_nbSamplesY; ++y)
		{
			double dx, dy;
			dx = double(x)/(m_nbSamplesX-1);
			dy = double(y)/(m_nbSamplesY-1);
			ImageType psd;
			psd.initItk(m_width, m_height, true);
			psd.for_all_pixels(m_function(dx, dy));
			m_explicitTPSD.addTexture(psd);
		}
	m_explicitTPSD.generate();
}

template<typename I>
const typename TPSD<I>::ImageType TPSD<I>::operator()(unsigned x, unsigned y) const
{
	double dx, dy;
	dx = double(x)/(m_nbSamplesX-1);
	dy = double(y)/(m_nbSamplesY-1);
	if(m_parallelogramCheat)
	{
		if(m_parallelogramOnX)
		{
			dx = dx-dy;
			if(dx<0)
				dx += 1;
		}
		else
		{
			dy = dy-dx;
			if(dy<0)
				dy += 1;
		}
	}
	ImageType psd;
	psd.initItk(m_width, m_height, true);
	psd.for_all_pixels(m_function(dx, dy));
	return psd;
}

}

}

#endif
