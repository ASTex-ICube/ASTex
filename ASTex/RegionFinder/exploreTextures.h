#ifndef __ASTEX_EXPLORETEXTURES__
#define __ASTEX_EXPLORETEXTURES__

#include "comparePixels.h"
#include <memory>

namespace ASTex
{

namespace RegionFinder
{

template<typename I, typename MASK_TYPE>
class RegionDistanceBase
{	
public:

	using RDOptionType = typename ComparePixelsBase<I, MASK_TYPE>::RDOptionBase;

	RegionDistanceBase():
		origin(),
		distance(0),
		rdOptions()
	{}

	RegionDistanceBase(const RegionDistanceBase& other)
		: origin(other.origin),
		  distance(other.distance),
		  rdOptions()
	{
		rdOptions = other.rdOptions;
	}

	virtual RegionDistanceBase& operator=(const RegionDistanceBase& other)
	{
		origin = other.origin;
		distance = other.distance;
		rdOptions = other.rdOptions;
		return *this;
	}

	virtual ~RegionDistanceBase()
	{}

	PixelPos origin;
	double distance;
	std::shared_ptr<RDOptionType> rdOptions; //maybe use std::unique_ptr instead?
};

template<typename I, typename MASK_TYPE>
class CompareRegionScoreBase
{	
public:

	using RegionDistanceType = RegionDistanceBase<I, MASK_TYPE>;

	CompareRegionScoreBase() {}
	bool operator()(const RegionDistanceType &object, const RegionDistanceType &other) const
	{
		return object.distance < other.distance;
	}
};

template<typename I, typename MASK_TYPE = float>
/**
 * @brief The ExploreTextureBase class serves as a base class for texture explorers.
 * An ExploreTexture class only needs to implement the operator() but may have other methods used for tweaks.
 */
class ExploreTextureBase
{
public:

	using RegionDistanceType = RegionDistanceBase<I, MASK_TYPE>;
	using CompareType = CompareRegionScoreBase<I, MASK_TYPE>;
	using PriorityQueueRFBase = std::priority_queue<RegionDistanceType, std::vector<RegionDistanceType>, CompareType>;
	using ImageMask = ImageGray<MASK_TYPE>;
	ExploreTextureBase() :
		m_forbiddenDistanceEpsilon(0)
	{}

	/**
	 * @brief setForbiddenScoreEpsilon is used to set a forbidden score under which a region is rejected.
	 * This function provides an easy way to detect and reject regions in the target
	 * that are identical or too close to the source.
	 * set epsilon to 0 or under to ignore this feature.
	 * @param epsilon distance under which a region is rejected.
	 */
	void setForbiddenDistanceEpsilon(double epsilon)
	{
		m_forbiddenDistanceEpsilon = epsilon;
	}

	/**
	 * @brief operator () iterates over targetTexture using compareRegionToTexture() (check examples).
	 *
	 * @return
	 */
	virtual PriorityQueueRFBase operator()(	const I& sourceRegion, const I& targetTexture,
												const ImageMask& sourceMask,
												const ImageMask& targetMask,
												ComparePixelsBase<I, MASK_TYPE>& comparePixels,
												unsigned nbRegionsKept, bool periodicity=false) const = 0;

protected:
	/**
	 * @brief compareRegionToTexture compares @param sourceRegion to @param targetTexture,
	 * in the position @param origin using @param comparePixels to compare pixels.
	 * @param sourceMask and @param targetMask define out of bounds regions respectively in the sourceRegion and in the
	 * targetTexture. @param periodicity defines if targetTexture and targtMask must be considered periodic or not.
	 * sourceRegion shall never be periodic.
	 * @return
	 */
	double compareRegionToTexture(	const I& sourceRegion, const I& targetTexture,
									const ImageMask& sourceMask, const ImageMask& targetMask,
									ComparePixelsBase<I, MASK_TYPE>& comparePixels,
									PixelPos origin, bool periodicity=false) const
	{
		assert(sourceRegion.is_initialized() && targetTexture.is_initialized()
			   && sourceMask.is_initialized() && targetMask.is_initialized() &&
			   "compareFullRegions: Some textures or masks are not initialized.");
		assert(sourceRegion.width()<targetTexture.width() && sourceRegion.height()<targetTexture.height() &&
			   "compareFullRegions: Source region must be smaller in size than the target texture.");
		assert(sourceRegion.size()==sourceMask.size() && targetTexture.size() == targetMask.size() &&
			   "compareFullRegions: Each mask needs to have the size of its corresponding texture.");
		bool suited = true;
		bool stopFlag = false;
		double distance = 0;
		PixelPos sourcePosition, targetPosition;
		for(unsigned x=0; x<unsigned(sourceRegion.width()) && suited && !stopFlag; ++x)
			for(unsigned y=0; y<unsigned(sourceRegion.height()) && suited && !stopFlag; ++y)
			{
				sourcePosition[0] = x;
				sourcePosition[1] = y;
				targetPosition[0] = periodicity ? (origin[0] + x)%targetTexture.width() : origin[0] + x;
				targetPosition[1] = periodicity ? (origin[1] + y)%targetTexture.height() : origin[1] + y;
				typename ImageMask::PixelType sourceMaskPix = sourceMask.pixelAbsolute(sourcePosition);
				typename ImageMask::PixelType targetMaskPix;
				//alpha map space: x, y. Original tile space: patchOrigin+(x, y). expTile space: origin+(x, y).
				if(sourceMaskPix>0)
				{
					if(!periodicity && (targetPosition[0] >= targetTexture.width()
										|| targetPosition[1] >= targetTexture.height()
										|| targetPosition[0] < 0
										|| targetPosition[1] < 0
						|| (targetMaskPix=targetMask.pixelAbsolute(targetPosition))==0) )
					{
						suited = false;
					}
					else
					{
						targetMaskPix=targetMask.pixelAbsolute(targetPosition);
						double comparasionResult = comparePixels(	sourceRegion,	targetTexture,
																	sourceMask,		targetMask,
																	sourcePosition,	targetPosition,
																	stopFlag,
																	periodicity);
						if(m_forbiddenDistanceEpsilon>0 && comparasionResult<m_forbiddenDistanceEpsilon)
							suited = false;
						else
							distance += comparasionResult * sourceMaskPix * targetMaskPix;
					}
				}
			}
		if(!suited)
			distance=std::numeric_limits<double>::infinity();
		return distance;
	}

	double			m_forbiddenDistanceEpsilon;
};

/**
 * @brief ExploreTextureExhaustive used to explore the entire target texture in order to find the very best matches.
 * Can be too slow for practical uses. Use ExploreTexturePyramid for a quicker exploration.
 */
template<typename I, typename MASK_TYPE = float>
class ExploreTextureExhaustive : public ExploreTextureBase<I, MASK_TYPE>
{
public:

	using PriorityQueueType = typename ExploreTextureBase<I, MASK_TYPE>::PriorityQueueRFBase;
	using RegionDistanceType = typename ExploreTextureBase<I, MASK_TYPE>::RegionDistanceType;
	using ImageMask = typename ExploreTextureBase<I, MASK_TYPE>::ImageMask;

	ExploreTextureExhaustive() :
		ExploreTextureBase<I, MASK_TYPE>()
	{}

	PriorityQueueType operator()(	const I& sourceRegion, const I& targetTexture,
									const ImageMask& sourceMask, const ImageMask& targetMask,
									ComparePixelsBase<I, MASK_TYPE>& comparePixels,
									unsigned nbRegionsKept, bool periodicity=false) const
	{
		PriorityQueueType pqc;
		double distance;
		PixelPos origin;
		for(origin[0]=0; origin[0]<targetTexture.width(); ++origin[0])
			for(origin[1]=0; origin[1]<targetTexture.height(); ++origin[1])
			{
				distance = ExploreTextureBase<I, MASK_TYPE>::compareRegionToTexture (	sourceRegion,	targetTexture,
																						sourceMask,		targetMask,
																						comparePixels,
																						origin, periodicity);
				RegionDistanceType rdt;
				rdt.origin = origin;
				rdt.distance = distance;
				rdt.rdOptions = comparePixels._rdOptions();
				pqc.push(rdt);
			}
		while(pqc.size()>nbRegionsKept)
			pqc.pop();
		return pqc;
	}
};

/**
 * @brief ExploreTextureSampler only explores regions randomly selected by a sampler.
 */
template<typename I, typename MASK_TYPE = float>
class ExploreTextureSampler : public ExploreTextureBase<I, MASK_TYPE>
{
public:

	using PriorityQueueType = typename ExploreTextureBase<I, MASK_TYPE>::PriorityQueueRFBase;
	using RegionDistanceType = typename ExploreTextureBase<I, MASK_TYPE>::RegionDistanceType;
	using ImageMask = typename ExploreTextureBase<I, MASK_TYPE>::ImageMask;

	ExploreTextureSampler(Stamping::SamplerBase &sampler) :
		ExploreTextureBase<I, MASK_TYPE>(),
		m_sampler(sampler)
	{}

	PriorityQueueType operator()(	const I& sourceRegion, const I& targetTexture,
									const ImageMask& sourceMask,
									const ImageMask& targetMask,
									ComparePixelsBase<I, MASK_TYPE>& comparePixels,
									unsigned nbRegionsKept, bool periodicity=false) const
	{
		PriorityQueueType pqc;
		double score;
		PixelPos origin;
		std::vector<Eigen::Vector2f> samples = m_sampler.generate();
		ImageGrayb hitmap;
		hitmap.initItk(targetMask.width(), targetMask.height(), true);
		for(auto sample : samples)
		{
			origin[0]=sample[0]*(targetTexture.width()-sourceRegion.width());
			origin[1]=sample[1]*(targetTexture.height()-sourceRegion.height());
			if(!hitmap.pixelAbsolute(origin[0], origin[1]))
			{
				hitmap.pixelAbsolute(origin[0], origin[1])=true;
				score = compareRegionToTexture (	sourceRegion,	targetTexture,
													sourceMask,		targetMask,
													comparePixels,
													origin, periodicity);
				RegionDistanceType rdt;
				rdt.origin = origin;
				rdt.distance = score;
				rdt.rdOptions = comparePixels._rdOptions();
				pqc.push(rdt);
			}
		}
		while(pqc.size()>nbRegionsKept)
			pqc.pop();
		return pqc;
	}

private:
	Stamping::SamplerBase &m_sampler;
};

/**
 * @brief ExploreTexturePyramid explores textures using a mipmap pyramid by first comparing the textures on
 * lower resolutions. Used to drastically increase the search without sacrifying too much quality.
 *
 */
template<typename I, typename MASK_TYPE = float>
class ExploreTexturePyramid : public ExploreTextureBase<I, MASK_TYPE>
{
public:
	using Base = ExploreTextureBase<I, MASK_TYPE>;
	using PriorityQueueType = typename Base::PriorityQueueRFBase;
	using RegionDistanceType = typename Base::RegionDistanceType;
	using RDOptionsType = typename ComparePixelsBase<I, MASK_TYPE>::RDOptionBase;
	using ImageMask = typename ExploreTextureBase<I, MASK_TYPE>::ImageMask;

	ExploreTexturePyramid() :
		m_minMapSize(64),
		m_trimTargetMask(true)
	{}

	/**
	 * @brief setMaxMipmapLevel sets the minimum size (maximum mipmap depth) of the target texture for the search.
	 */
	void setMinimumTargetMipmapSize(unsigned size)
	{
		m_minMapSize = size;
	}

	/**
	 * @brief setTrimTargetMask sets if the target Mask should be trimmed on each resolution, or expanded.
	 * true by default, because it causes less issues with avoided regions suddenly becoming unavoided on lower res.
	 * The problem comes from the mipmapping of the masks causing regions to potentially bleed into each others.
	 */
	void setTrimTargetMask(bool trim)
	{
		m_trimTargetMask = trim;
	}

	PriorityQueueType operator()(	const I& sourceRegion, const I& targetTexture,
									const ImageMask& sourceMask,
									const ImageMask& targetMask,
									ComparePixelsBase<I, MASK_TYPE>& comparePixels,
									unsigned nbRegionsKept, bool periodicity=false) const
	{
		//Precondition for this operator: each of those arguments must be their corresponding mipmap's first map.
		//If possible, an assert should verify this.
		double infinity = std::numeric_limits<double>::infinity();
		Mipmap<I> sourceRegionMipmap(sourceRegion), targetTextureMipmap(targetTexture);
		Mipmap<ImageMask> sourceMaskMipmap(sourceMask), targetMaskMipmap(targetMask);
		sourceRegionMipmap.generate();
		targetTextureMipmap.generate();
		sourceMaskMipmap.generate();
		targetMaskMipmap.generate();

		for(unsigned k=0; k<targetTextureMipmap.numberMipmapsWidth(); ++k)
			targetMaskMipmap.mipmap(k, k).for_all_pixels([&] (typename ImageMask::PixelType &pix)
			{
				if(pix < 1.0)
					pix=m_trimTargetMask ? 0.0 : 1.0;
			});

		PriorityQueueType regionComparatorPQ;

		int maxK = int(std::max(double(0), targetTextureMipmap.numberMipmapsWidth()-std::log(m_minMapSize)));
		for(int k=maxK; k>=0; --k) //^ back is where the tile with identity transform is (rember if modified)
		{
			PriorityQueueType ccPQCopy;
			if(k!=maxK)
			{
				while(!regionComparatorPQ.empty())
				{
					const I& lowResSourceRegion = sourceRegionMipmap.mipmap(k, k);
					RegionDistanceType rdtNew;
					const RegionDistanceType &rdtTop = regionComparatorPQ.top();
					PixelPos newOrigin = rdtTop.origin;
					double distances[4], bestDistance=infinity;
					std::shared_ptr<RDOptionsType> bestRdOptions = nullptr;
					unsigned bestIndex=0;
					PixelPos origins[4];
					newOrigin[0] *= 2;
					newOrigin[1] *= 2;
					origins[0] = newOrigin;
					++newOrigin[0];
					origins[1] = newOrigin;
					++newOrigin[1];
					origins[2] = newOrigin;
					--newOrigin[0];
					origins[3] = newOrigin;
					for(unsigned i=0; i<4; ++i)
					{
						distances[i] = ExploreTextureBase<I, MASK_TYPE>::compareRegionToTexture(lowResSourceRegion,
																							 targetTextureMipmap.mipmap(k, k),
																							 sourceMaskMipmap.mipmap(k, k),
																							 targetMaskMipmap.mipmap(k, k),
																							 comparePixels, origins[i], periodicity);
						if(distances[i] < bestDistance)
						{
							bestIndex = i;
							bestDistance = distances[i];
							bestRdOptions = comparePixels._rdOptions();
						}
					}
					if(bestDistance == infinity)
					{	//v the process may land out of bounds if the level above (resolution bellow) has an odd size
						//compared to the previously treated resolution, even though this is kind of unusual
						--newOrigin[1];
						double finalDistances[2] = {infinity, infinity};
						PixelPos finalOrigins[2] = {newOrigin, newOrigin};
						if(	sourceRegionMipmap.mipmap(k, k).width()%2 != 0 &&
							sourceRegionMipmap.mipmap(k, k).height()%2 == 0)
						{
							--finalOrigins[0][0];
							--finalOrigins[1][0];
							++finalOrigins[1][1];
						}
						else if(sourceRegionMipmap.mipmap(k, k).width()%2 == 0 &&
								sourceRegionMipmap.mipmap(k, k).height()%2 != 0)
						{
							--finalOrigins[0][1];
							--finalOrigins[1][1];
							++finalOrigins[1][0];
						}
						if(sourceRegionMipmap.mipmap(k, k).width()%2 != 0 &&
								sourceRegionMipmap.mipmap(k, k).height()%2 != 0)
						{
							--finalOrigins[0][0];
							--finalOrigins[0][1];
							finalDistances[0] = ExploreTextureBase<I, MASK_TYPE>::compareRegionToTexture(sourceRegionMipmap.mipmap(k, k),
																										targetTextureMipmap.mipmap(k, k),
																										sourceMaskMipmap.mipmap(k, k),
																										targetMaskMipmap.mipmap(k, k),
																										comparePixels, finalOrigins[0],
																										periodicity);
						}
						else
						{
							for(unsigned i=0; i<2; ++i)
							{
								finalDistances[i] = ExploreTextureBase<I, MASK_TYPE>::compareRegionToTexture(sourceRegionMipmap.mipmap(k, k),
																											targetTextureMipmap.mipmap(k, k),
																											sourceMaskMipmap.mipmap(k, k),
																											targetMaskMipmap.mipmap(k, k),
																											comparePixels, finalOrigins[i],
																											 periodicity);
								if(finalDistances[i] < bestDistance)
								{
									bestIndex = i;
									origins[i] = finalOrigins[i];
									distances[i] = finalDistances[i];
									bestDistance = distances[i];
									bestRdOptions = comparePixels._rdOptions();
								}
							}
						}
					}
					regionComparatorPQ.pop();
					rdtNew.distance = bestDistance;
					rdtNew.origin = origins[bestIndex];
					if(bestRdOptions != nullptr)
					{
						rdtNew.rdOptions = bestRdOptions;
					}
					else
						rdtNew.rdOptions = nullptr;
					ccPQCopy.push(rdtNew);
				}
				regionComparatorPQ = ccPQCopy;
			}
			else
			{	//first step: a proportional number of contents per texture.
				ExploreTextureExhaustive<I, MASK_TYPE> exploreTextureExhaustive;
				regionComparatorPQ = exploreTextureExhaustive(	sourceRegionMipmap.mipmap(k, k), targetTextureMipmap.mipmap(k, k),
																sourceMaskMipmap.mipmap(k, k), targetMaskMipmap.mipmap(k, k),
																comparePixels, nbRegionsKept, periodicity);
			}
			while(regionComparatorPQ.size() > nbRegionsKept)
				regionComparatorPQ.pop();
		}
		return regionComparatorPQ;
	}

private:
	unsigned		m_minMapSize;
	bool			m_trimTargetMask;
};

}

}

#endif
