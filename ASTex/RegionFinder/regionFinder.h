#ifndef __ASTEX_REGIONFINDER__
#define __ASTEX_REGIONFINDER__

#include "exploreTextures.h"

//Pixel comparators

namespace ASTex
{

namespace RegionFinder
{

/**
 * @brief findBestRegions Function used whenever matches have to be found between a source and a target.
 * The source is always smaller than the target and is the region the function will try to find a match for.
 * Masks can be used to discriminate the search. If you do not care about that, provide masks that are 1-valued
 * and of the size of their owner.
 * @param sourceRegion the source region to find matches for.
 * @param targetTexture the target texture to find matches in.
 * @param sourceMask The mask of the source, used to ignore texels or weight them down during the search.
 * For each texel:	1 means the texel is not to be ignored.
					Between 0 and 1 excluded means the texel is weighed by the same number.
					0 means the texel is to be ignored, and the same texel found in the target will always be valid.
					Above 1 and bellow 0 is undefined behaviour.
 * @param targetMask The mask of the target, used to avoid texels during the search.
 * For each texel:	1 means the texel is not to be avoided.
 *					Between 0 and 1 excluded means the texel is weighed by the same number (useless in most cases).
 *					0 means the texel is to be avoided and will be considered out of bounds.
 *					Above 1 and bellow 0 is undefined behaviour.
 * @param exploreTexture Function defining how the target will be explored. Typical uses are already coded.
 * @param comparePixels Function defining how the source will be compared to the target. Works per-pixel or per-region.
 * @param nbRegionsKept Number of regions to keep in the priority queue returned.
 * @param periodicity If the target is considered to be a periodic image (a.k.a repeating on each side) or not.
 * @return
 */
template<typename I, typename MASK_TYPE=float>
typename ExploreTextureBase<I, MASK_TYPE>::PriorityQueueRFBase
						findBestRegions(const I& sourceRegion, const I& targetTexture,
										const typename ExploreTextureBase<I, MASK_TYPE>::ImageMask& sourceMask,
										const typename ExploreTextureBase<I, MASK_TYPE>::ImageMask& targetMask,
										const ExploreTextureBase<I, MASK_TYPE>& exploreTexture,
										ComparePixelsBase<I, MASK_TYPE>& comparePixels,
										unsigned nbRegionsKept, bool periodicity=false)
{
	assert(sourceRegion.is_initialized() && targetTexture.is_initialized()
		   && sourceMask.is_initialized() && targetMask.is_initialized() &&
		   "findBestRegion: Some textures or masks are not initialized.");
	assert(sourceRegion.width()<targetTexture.width() && sourceRegion.height()<targetTexture.height() &&
		   "findBestRegion: Source region must be smaller in size than the target texture.");
	assert(sourceRegion.size()==sourceMask.size() && targetTexture.size() == targetMask.size() &&
		   "findBestRegion: Each mask needs to have the size of its corresponding texture.");

	return exploreTexture(sourceRegion, targetTexture, sourceMask, targetMask, comparePixels, nbRegionsKept, periodicity);
}

}

}

#endif
