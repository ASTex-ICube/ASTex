#ifndef __ASTEX_COMPAREPIXELS__
#define __ASTEX_COMPAREPIXELS__

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <cstring>

#include "ASTex/image_gray.h"
#include "ASTex/image_rgb.h"
#include "ASTex/image_rgba.h"
#include "ASTex/utils.h"
#include "ASTex/Stamping/sampler.h"
#include "ASTex/mipmap.h"

namespace ASTex
{

namespace RegionFinder
{

using PixelPos = itk::Index<2>;

/**
 * @brief ComparePixelsBase base class for ComparePixels classes.
 * I is the texture type and MASK_TYPE the sub-type of the mask.
 */
template<typename I, typename MASK_TYPE=float>
class ComparePixelsBase
{
public:

	class RDOptionBase
	{
	public:

		RDOptionBase()
		{}

		virtual ~RDOptionBase()
		{}

		virtual RDOptionBase& operator=(const RDOptionBase& other)
		{
			(void) other;
			return *this;
		}
	};

	using ImageMask = ImageGray<MASK_TYPE>;

	ComparePixelsBase() :
		m_rdOptions(nullptr)
	{}

	virtual ~ComparePixelsBase()
	{}

	virtual double operator()(	const I& sourceRegion, const I& targetTexture,
								const ImageMask& sourceMask,
								const ImageMask& targetMask,
								PixelPos sourcePosition, PixelPos targetPosition,
								bool &stopFlag,
								bool periodicity=false) = 0;

	/**
	 * @brief _rdOptions is only used to transfer options related to the result of a search to the ExploreTexture class.
	 */
	std::shared_ptr<RDOptionBase> _rdOptions() {return m_rdOptions;}

protected:
	std::shared_ptr<RDOptionBase> m_rdOptions;
};

/**
 * @brief ComparePixelsMSE Compares regions by their mean squared error.
 * A good toy example to make your own ComparePixels class.
 */
template<typename I, typename MASK_TYPE=float>
class ComparePixelsMSE : public ComparePixelsBase<I, MASK_TYPE>
{
public:

	ComparePixelsMSE() :
		ComparePixelsBase<I, MASK_TYPE>()
	{}

	using ImageMask = typename ComparePixelsBase<I, MASK_TYPE>::ImageMask;

	double operator()(	const I& sourceRegion, const I& targetTexture,
						const ImageMask& sourceMask,
						const ImageMask& targetMask,
						PixelPos sourcePosition, PixelPos targetPosition,
						bool &stopFlag,
						bool periodicity=false)
	{
		(void) sourceMask;
		(void) targetMask;
		(void) stopFlag;
		return mse(sourceRegion, targetTexture,
				sourcePosition[0], sourcePosition[1],
				targetPosition[0], targetPosition[1], 0, periodicity);

	}
};

/**
 * @brief ComparePixelsMean Compares regions by their average values.
 * The closer, the smaller the distance (norm 2 for multi-channels).
 */
template<typename I, typename MASK_TYPE=float>
class ComparePixelsMean : public ComparePixelsBase<I, MASK_TYPE>
{
public:

	ComparePixelsMean() :
		ComparePixelsBase<I, MASK_TYPE>()
	{}

	using ImageMask = typename ComparePixelsBase<I, MASK_TYPE>::ImageMask;

	double operator()(	const I& sourceRegion, const I& targetTexture,
						const ImageMask& sourceMask,
						const ImageMask& targetMask,
						PixelPos sourceOrigin, PixelPos targetOrigin,
						bool &stopFlag,
						bool periodicity=false)
	{
		(void) sourceOrigin;
		(void) periodicity;

		double distance = 0;
		if(!stopFlag)
		{
			using PixelType = typename I::PixelType;
			using DataType = typename I::DataType;
			unsigned pixelSize = sizeof(PixelType) / sizeof(DataType);
			typename I::PixelType meanSource = I::zero(), meanTarget = I::zero();
			typename ImageMask::PixelType smTotal = ImageMask::zero(), tmTotal = ImageMask::zero();
			bool outOfBounds = false;
			for(unsigned x=0; x<unsigned(sourceRegion.width()) && !outOfBounds; ++x)
				for(unsigned y=0; y<unsigned(sourceRegion.height()) && !outOfBounds; ++y)
				{
					typename ImageMask::PixelType smPix = sourceMask.pixelAbsolute(x, y);
					if(smPix>0)
					{
						PixelPos targetCoordinates;
						if(periodicity)
						{
							targetCoordinates[0] = (targetOrigin[0]+x)%targetTexture.width();
							targetCoordinates[1] = (targetOrigin[1]+y)%targetTexture.height();
						}
						else
						{
							targetCoordinates[0] = targetOrigin[0]+x;
							targetCoordinates[1] = targetOrigin[1]+y;
						}
						typename ImageMask::PixelType tmPix;
						if((targetCoordinates[0]<0 || targetCoordinates[1]<0 ||
							targetCoordinates[0]>=targetTexture.width() ||
							targetCoordinates[1]>=targetTexture.height() ||
							(tmPix=targetMask.pixelAbsolute(targetCoordinates))==0) )
							outOfBounds=true;
						else
						{
							tmTotal += tmPix;
							smTotal += smPix;
							meanSource += sourceRegion.pixelAbsolute(x, y)*smPix;
							meanTarget += targetTexture.pixelAbsolute(x, y)*tmPix;
						}
					}
				}
			stopFlag = true;
			if(!outOfBounds)
			{
				for(unsigned x=0; x<unsigned(sourceRegion.width()); ++x)
					for(unsigned y=0; y<unsigned(sourceRegion.height()); ++y)
					{
						typename ImageMask::PixelType smPix = sourceMask.pixelAbsolute(x, y);
						if(smPix>0)
						{
							PixelPos targetCoordinates;
							typename ImageMask::PixelType tmPix;
							targetCoordinates[0] = targetOrigin[0]+x;
							targetCoordinates[1] = targetOrigin[1]+y;
							assert(targetMask.pixelAbsolute(targetCoordinates)==0);
						}
					}
				meanSource = meanSource * (1.0/smTotal);
				meanTarget = meanTarget * (1.0/tmTotal);
				DataType *meanSourceDataType = reinterpret_cast<DataType *>(&meanSource);
				DataType *meanTargetDataType = reinterpret_cast<DataType *>(&meanTarget);
				for(unsigned i=0; i<pixelSize; ++i)
				{
					distance +=		(meanSourceDataType[i] - meanTargetDataType[i])
								*	(meanSourceDataType[i] - meanTargetDataType[i]);
				}
			}
			else
				distance = std::numeric_limits<double>::infinity();
		}
		return distance;
	}
};

/**
 * @brief ComparePixelsScalarProduct Compares regions by projecting one onto another, while normalizing these regions.
 * Note that to get the norms after the search, the RegionDistanceOption (RdOption) has to be read.
 * Also note that the distance is a negative number; the largest projection is -(number of channels of the texture type)
 */
template<typename I, typename MASK_TYPE=float>
class ComparePixelsScalarProduct : public ComparePixelsBase<I, MASK_TYPE>
{
public:

	using RDOptionBase = typename ComparePixelsBase<I, MASK_TYPE>::RDOptionBase;
	using PixelType = typename I::PixelType;
	using DataType = typename I::DataType;
	using MaskPixelType = MASK_TYPE;
	using ImageMask = typename ComparePixelsBase<I, MASK_TYPE>::ImageMask;

	class RDOptionScalarProduct : public RDOptionBase
	{
	public:
		RDOptionScalarProduct() :
			RDOptionBase(),
			sourceNorm(),
			targetNorm()
		{}

		virtual RDOptionScalarProduct& operator=(const RDOptionScalarProduct& other)
		{
			sourceNorm = other.sourceNorm;
			targetNorm = other.targetNorm;
			return *this;
		}

		virtual ~RDOptionScalarProduct()
		{}

		PixelType sourceNorm;
		PixelType targetNorm;
	};

	ComparePixelsScalarProduct() :
		ComparePixelsBase<I, MASK_TYPE>()
	{}

	double operator()(	const I& sourceRegion, const I& targetTexture,
						const ImageMask& sourceMask,
						const ImageMask& targetMask,
						PixelPos sourceOrigin, PixelPos targetOrigin,
						bool &stopFlag,
						bool periodicity=false)
	{
		(void) sourceOrigin;
		(void) periodicity;

		double distance = 0;
		if(!stopFlag)
		{
			unsigned pixelSize = sizeof(PixelType) / sizeof(DataType);
			bool outOfBounds = false;
			std::shared_ptr<RDOptionScalarProduct> options = std::make_shared<RDOptionScalarProduct>();
			this->m_rdOptions = options;
			PixelType sourceNorm = I::zero(), targetNorm = I::zero();
			DataType *sourceNormData = reinterpret_cast<DataType *>(&sourceNorm);
			DataType *targetNormData = reinterpret_cast<DataType *>(&targetNorm);
			for(unsigned x=0; x<unsigned(sourceRegion.width()) && !outOfBounds; ++x)
				for(unsigned y=0; y<unsigned(sourceRegion.height()) && !outOfBounds; ++y)
				{
					MaskPixelType smPix = sourceMask.pixelAbsolute(x, y);
					if(smPix>0)
					{
						PixelPos targetCoordinates;
						if(periodicity)
						{
							targetCoordinates[0] = (targetOrigin[0]+x)%targetTexture.width();
							targetCoordinates[1] = (targetOrigin[1]+y)%targetTexture.height();
						}
						else
						{
							targetCoordinates[0] = targetOrigin[0]+x;
							targetCoordinates[1] = targetOrigin[1]+y;
						}
						typename ImageMask::PixelType tmPix;
						if((targetCoordinates[0]<0 || targetCoordinates[1]<0 ||
							targetCoordinates[0]>=targetTexture.width() ||
							targetCoordinates[1]>=targetTexture.height() ||
							(tmPix=targetMask.pixelAbsolute(targetCoordinates))==0) )
							outOfBounds=true;
						else
						{
							PixelType targetPix = targetTexture.pixelAbsolute(targetCoordinates);
							PixelType sourcePix = sourceRegion.pixelAbsolute(x, y);
							for(unsigned i=0; i<pixelSize; ++i)
							{
								sourceNormData[i] += sourcePix[i] * sourcePix[i] * smPix * tmPix;
								targetNormData[i] += targetPix[i] * targetPix[i] * smPix * tmPix;
							}
						}
					}
				}
			stopFlag = true;
			if(!outOfBounds)
			{
				for(unsigned i=0; i<pixelSize; ++i)
				{
					sourceNormData[i] = std::sqrt(double(sourceNormData[i]));
					targetNormData[i] = std::sqrt(double(targetNormData[i]));
				}
				options->sourceNorm = sourceNorm;
				options->targetNorm = targetNorm;
				for(unsigned x=0; x<unsigned(sourceRegion.width()) && !outOfBounds; ++x)
					for(unsigned y=0; y<unsigned(sourceRegion.height()) && !outOfBounds; ++y)
					{
						MaskPixelType smPix = sourceMask.pixelAbsolute(x, y);
						if(smPix>0)
						{
							PixelPos targetCoordinates;
							targetCoordinates[0] = periodicity ? (targetOrigin[0]+x)%targetTexture.width() : targetOrigin[0]+x;
							targetCoordinates[1] = periodicity ? (targetOrigin[1]+y)%targetTexture.height() : targetOrigin[1]+y;
							MaskPixelType tmPix = targetMask.pixelAbsolute(targetCoordinates);
							PixelType targetPix = targetTexture.pixelAbsolute(targetCoordinates);
							PixelType sourcePix = sourceRegion.pixelAbsolute(x, y);
							DataType *sourcePixData = reinterpret_cast<DataType *>(&sourcePix);
							DataType *targetPixData = reinterpret_cast<DataType *>(&targetPix);
							double distancePix = 0;
							for(unsigned i=0; i<pixelSize; ++i)
							{
								double distanceChannel = 0;
								distanceChannel =	sourcePixData[i] * smPix * tmPix * (1.0/sourceNormData[i]) *
													targetPixData[i] * smPix * tmPix * (1.0/targetNormData[i]);
//								distanceChannel *= distanceChannel;
								distancePix += distanceChannel;
							}
							distance -= distancePix;
						}
					}
			}
			else
			{
				distance = std::numeric_limits<double>::infinity();

			}
		}
		return distance;
	}
};

}

}

#endif
