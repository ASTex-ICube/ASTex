#ifndef __CTEXCH_MULTI_TEXTURE_PROCESSOR_H__
#define __CTEXCH_MULTI_TEXTURE_PROCESSOR_H__

#include <random>
#include <algorithm>

#include "patchProcessor.h"
#include "ASTex/Stamping/sampler.h"
#include "ASTex/RegionFinder/regionFinder.h"
#include "ASTex/texture_pool.h"



namespace ASTex
{

namespace ContentExchange
{

template<typename I>
/**
 * @brief The MultiTextureProcessor class is a mod made for
 * processing and storing multi-texture content exchange elements.
 * Feed it a tile, 0 to n textures (precise which are periodic textures)
 * and 0 to m transformations (Eigen matrices),
 * choose how you want patches and contents to be computed, then use
 * generate() to get an output tiling, or export the class for online rendering.
 * This class has a gradient mode, for correctly handling rotations in the gradient domain.
 */
class MultiTextureProcessor : public ASTex::ContentExchange::PatchProcessor<I>
{
public:

	using ImageType = I;
	using PixelType = typename ImageType::PixelType;
	using DataType = typename ImageType::DataType;
	using PriorityQueueType = typename RegionFinder::ExploreTextureBase<I, double>::PriorityQueueRFBase;
	using RegionDistanceType = RegionFinder::RegionDistanceBase<I, double>;

	MultiTextureProcessor();
	MultiTextureProcessor(const ImageType& image);

	~MultiTextureProcessor();

	void setTexturePool(const TexturePool<ImageType, double> texturePool)
	{
		m_transformedTexturePool = texturePool;
	}

	void contents_explorePools();

	Mipmap<I> generateFromExemplar(const I& exemplar);

private:

	TexturePool<ImageType, double> m_transformedTexturePool;
};

template<typename I>
MultiTextureProcessor<I>::MultiTextureProcessor():
	PatchProcessor<ImageType>()
{}

template<typename I>
MultiTextureProcessor<I>::MultiTextureProcessor(const ImageType& tile):
	PatchProcessor<ImageType>(tile)
{}

template<typename I>
MultiTextureProcessor<I>::~MultiTextureProcessor()
{}

template<typename I>
void MultiTextureProcessor<I>::contents_explorePools()
{
	assert(this->m_patches.size()>0 &&
		   "MultiTextureProcessor::contents_explorePools: a patch initializer must be called before being able to choose contents");
	assert(this->m_mipmapMode != NO_FILTER &&
			"MultiTextureProcessor::contents_explorePools: Mipmap mode cannot be NO_FILTER even when ignoring mip-maps for this function");
	assert(this->m_patches[0].nbContents()>0 &&
			"MultiTextureProcessor::contents_explorePools: default contents must have been set.");
	srand(this->m_seed);
	Mipmap<ImageType> inputTileMipmap(this->m_tile);
	inputTileMipmap.generate();

	PriorityQueueType contentComparatorPQ;
	using PairDistanceIndexType = std::pair<RegionDistanceType, unsigned>;

	for(unsigned p=0; p<this->m_patches.size(); ++p)
	{	
		std::vector<PairDistanceIndexType> texContentComparatorVector;
		//We have to do some cleansing before we're able to use the data computed in the content exchange process.

		Patch<ImageType> &patch=this->m_patches[p];
		const Content<ImageType> &defaultContent = patch.contentAt(0);
		Mipmap<ImageAlphad> cleansedPatchMipmap;
		ImageAlphad cleansedPatch;
		cleansedPatch.initItk(this->m_tile.width(), this->m_tile.height(), true);
		//^ v (kind of silly, but this is a way to get the right number of reductions from Mipmap)
		cleansedPatchMipmap.setTexture(cleansedPatch);
		cleansedPatchMipmap.generate();

		I cleansedContent;
		cleansedContent.initItk(this->m_tile.width(), this->m_tile.height());
		cleansedContent.copy_pixels(this->m_tile);

		// v now we treat the patch mip-map that we allocated. We need to ignore texels in the center of the patch.
		cleansedPatch.initItk(patch.alphaMap().width(), patch.alphaMap().height());
		cleansedPatch.for_all_pixels([&] (ImageAlphad::PixelType &pix, int x, int y)
		{
			ImageAlphad::PixelType pixPatch = patch.alphaMap().pixelAbsolute(x, y);
			pix = pixPatch >= 1 ? 0 : pixPatch; //We keep only texels which are close to the border on each resolution.
			//This is not very portable because it requires us to always call patches_dilate() before using this function
			//But it is just so convenient and elegant
		});

		// v finally, we treat the default content. The problem with the default content mip-map is that texels are
		//blackened near the border, which adds a huge bias to the region search.

		const I& defaultContentMap = defaultContent.texture();
		cleansedContent.initItk(defaultContentMap.width(), defaultContentMap.height(), true);
		cleansedContent.for_all_pixels([&] (PixelType &pix, int x, int y)
		{
			ImageAlphad::PixelType pixPatch = patch.alphaMap().pixelAbsolute(x, y);
			if(pixPatch>0)
			{
				pix = defaultContentMap.pixelAbsolute(x, y) * (1.0/pixPatch);
			}
		});
		unsigned nbKeptContentsPerTexture = this->nbContentsPerPatch()*2; //< I don't know what to put here

		for(unsigned i=0; i<m_transformedTexturePool.size(); ++i)
		{
			const typename TexturePool<ImageType, double>::TransformedTexture &transformedTexture = m_transformedTexturePool[i];
			RegionFinder::ExploreTextureExhaustive<ImageType, double> exploreTexture;
			//exploreTexture.setForbiddenDistanceEpsilon(0.000001);
			RegionFinder::ComparePixelsMSE<ImageType, double> comparePixels;
			contentComparatorPQ = RegionFinder::findBestRegions(	cleansedContent, transformedTexture.texture,
																	cleansedPatchMipmap.texture(), transformedTexture.boundaries,
																	exploreTexture, comparePixels,
																	nbKeptContentsPerTexture, transformedTexture.periodicity);
			while(contentComparatorPQ.size() > 0)
			{
				RegionDistanceType rdt = contentComparatorPQ.top();
				contentComparatorPQ.pop();
				texContentComparatorVector.push_back(std::make_pair(rdt, i));
			}
		}
		std::sort(texContentComparatorVector.begin(), texContentComparatorVector.end(),
				  [](const PairDistanceIndexType& object, const PairDistanceIndexType &other)
		{
			return object.first.distance < other.first.distance;
		});
		//Now we can build the contents out of the remains!
		for(unsigned i=0; i<this->nbContentsPerPatch(); ++i)
		{
			RegionDistanceType contentComparator = texContentComparatorVector[i].first;
			I contentTile;
			const I& selectedTexture = m_transformedTexturePool[texContentComparatorVector[i].second].texture;
			contentTile.initItk(this->m_tile.width(), this->m_tile.height(), true);
			for(int x=0; x<patch.alphaMap().width(); ++x)
				for(int y=0; y<patch.alphaMap().height(); ++y)
				{
					unsigned xPatch = (x+patch.originAt(0, 0)[0])%this->m_tile.width();
					unsigned yPatch = (y+patch.originAt(0, 0)[1])%this->m_tile.height();
					unsigned xContent = m_transformedTexturePool[texContentComparatorVector[i].second].periodicity ?
								(x+contentComparator.origin[0]) % selectedTexture.width()
						:		(x+contentComparator.origin[0]);
					unsigned yContent = m_transformedTexturePool[texContentComparatorVector[i].second].periodicity ?
								(y+contentComparator.origin[1]) % selectedTexture.height()
						:		(y+contentComparator.origin[1]);
					if(patch.alphaMap().pixelAbsolute(x, y)>0) //very important because some pixels may be out of bounds
						contentTile.pixelAbsolute(xPatch, yPatch) = selectedTexture.pixelAbsolute(xContent, yContent);
				}
			Content<ImageType> content(contentTile, patch);
			patch.addContent(content);
		}
	}
}

template<typename I>
Mipmap<I> MultiTextureProcessor<I>::generateFromExemplar(const I &exemplar)
{
	assert(this->m_patches.size()>0 && "PatchProcessor::generate: a patch initializer must be called before using generate()");
	ImageType output;
	output.initItk(this->m_outputWidth, this->m_outputHeight, true);

	ImageType scaledExemplar;
	scaledExemplar.initItk(this->m_outputWidth, this->m_outputHeight, true);
	scaledExemplar.for_all_pixels([&] (PixelType &pix, int x, int y)
	{
		double xd, yd;
		xd = double(x)/(scaledExemplar.width()-1) * (exemplar.width()-1);
		yd = double(y)/(scaledExemplar.height()-1) * (exemplar.height()-1);
		pix = bilinear_interpolation(exemplar, xd, yd, false);
	});

	ImageRGBd uvMap;
	uvMap.initItk(this->m_outputWidth, this->m_outputHeight, true);

	unsigned pixelSize = sizeof(PixelType)/sizeof(DataType);
	srand(this->m_seed);

	for(unsigned p=0; p<this->m_patches.size(); ++p)
	{
		Patch<ImageType> &patch = this->m_patches[p];
		PixelPos patchOrigin = patch.originAt(0, 0);
		PixelPos outputOrigin = patchOrigin;
		for(outputOrigin[0]=patchOrigin[0]-this->m_tile.width();
			outputOrigin[0]<unsigned(this->m_outputWidth);
			outputOrigin[0]+=this->m_tile.width())
			for(outputOrigin[1]=patchOrigin[1]-this->m_tile.height();
				outputOrigin[1]<unsigned(this->m_outputHeight);
				outputOrigin[1]+=this->m_tile.height())
			{
				bool entirelyOutOfBounds = true;
				//compute a content for this patch such that its region corresponds to
				//the region we try to find a best match for.
				//The reason for using a different alpha map than that of the patch is because patches may be cut
				//on the borders of the output.
				ImageAlphad emulatedAlpha;
				ImageType emulatedTile;
				emulatedAlpha.initItk(this->m_tile.width(), this->m_tile.height(), true);
				emulatedTile.initItk(this->m_tile.width(), this->m_tile.height(), true);
				for(int x=0; x<patch.alphaMap().width(); ++x)
					for(int y=0; y<patch.alphaMap().height(); ++y)
					{
						if(	outputOrigin[0]+x>=0 && outputOrigin[1]+y>=0
							&& outputOrigin[0]+x<scaledExemplar.width() && outputOrigin[1]+y<scaledExemplar.height()
							&& patch.alphaMap().pixelAbsolute(x, y)>0)
						{
							PixelPos tileCoordinates;
							PixelPos patchCoordinates;
							tileCoordinates[0] = (patchOrigin[0]+x)%emulatedTile.width();
							tileCoordinates[1] = (patchOrigin[1]+y)%emulatedTile.height();
							patchCoordinates[0] = outputOrigin[0]+x;
							patchCoordinates[1] = outputOrigin[1]+y;
							emulatedAlpha.pixelAbsolute(tileCoordinates)= patch.alphaMap().pixelAbsolute(x, y);
							emulatedTile.pixelAbsolute(tileCoordinates)= scaledExemplar.pixelAbsolute(patchCoordinates);
							entirelyOutOfBounds = false;
						}
					}
				if(!entirelyOutOfBounds)
				{
					using PairDistanceIndexType = std::pair<RegionDistanceType, unsigned>;
					std::vector<PairDistanceIndexType> texContentComparatorVector;
					Patch<ImageType> emulatedPatch; //creating a patch also automatically cuts empty pixels out
					emulatedPatch.setAlphaMap(emulatedAlpha, ISOTROPIC);
					Content<ImageType> emulatedContent(emulatedTile, emulatedPatch);
					for(unsigned i=0; i<m_transformedTexturePool.size(); ++i)
					{
						PriorityQueueType contentComparatorPQ;
						const typename TexturePool<ImageType, double>::TransformedTexture &transformedTexture = m_transformedTexturePool[i];
						PriorityQueueType subContentComparatorPQ;
						RegionFinder::ExploreTexturePyramid<ImageType, double> exploreTexture;
						exploreTexture.setTrimTargetMask(true);
						RegionFinder::ComparePixelsScalarProduct<ImageType, double> comparePixels;
						contentComparatorPQ = RegionFinder::findBestRegions(emulatedContent.texture(),
																			transformedTexture.texture,
																			emulatedPatch.alphaMap(),
																			transformedTexture.boundaries,
																			exploreTexture, comparePixels,
																			this->m_nbContentsPerPatch,
																			transformedTexture.periodicity);
						while(contentComparatorPQ.size() > 0)
						{
							RegionDistanceType rdt = contentComparatorPQ.top();
							contentComparatorPQ.pop();
							texContentComparatorVector.push_back(std::make_pair(rdt, i));
						}
					}
					std::sort(texContentComparatorVector.begin(), texContentComparatorVector.end(),
							  [](const PairDistanceIndexType& object, const PairDistanceIndexType &other)
					{
						return object.first.distance < other.first.distance;
					});
					//std::random_shuffle(texContentComparatorVector.begin(), texContentComparatorVector.begin()
					//					+std::min(texContentComparatorVector.size(), size_t(this->m_nbContentsPerPatch))
					//					);
					std::mt19937 rd_engine(1234567);
					std::shuffle(texContentComparatorVector.begin(), texContentComparatorVector.begin()
						+ std::min(texContentComparatorVector.size(), size_t(this->m_nbContentsPerPatch)),rd_engine);

					
					const PairDistanceIndexType& tex = texContentComparatorVector[0];
					using RDOptionsType = typename RegionFinder::ComparePixelsScalarProduct<ImageType, double>::RDOptionScalarProduct;
					std::shared_ptr<RDOptionsType> rdosp = std::dynamic_pointer_cast<RDOptionsType>(tex.first.rdOptions); //<here
					std::cout << "Best - origin: " << tex.first.origin << ", distance: " << tex.first.distance << ", texture pool id: " << tex.second << std::endl;
					I foundExpTexture;
					foundExpTexture.initItk(this->m_tile.width(), this->m_tile.height(), true);
					PixelPos emulatedPatchOrigin = emulatedPatch.originAt(0, 0);
					PixelPos bestTextureOrigin = texContentComparatorVector[0].first.origin;
					PixelPos foundTextureCoordinates;
					const typename TexturePool<I, double>::TransformedTexture &transformedTexture
							= m_transformedTexturePool[texContentComparatorVector[0].second];
					for(int x=0; x<emulatedPatch.alphaMap().width(); ++x)
						for(int y=0; y<emulatedPatch.alphaMap().height(); ++y)
						{
							if(emulatedPatch.alphaMap().pixelAbsolute(x, y)>0)
							{
								bool periodicity = transformedTexture.periodicity;
								PixelPos tileCoordinates;
								tileCoordinates[0] = (emulatedPatchOrigin[0]+x)%foundExpTexture.width();
								tileCoordinates[1] = (emulatedPatchOrigin[1]+y)%foundExpTexture.height();
								const I &bestTexture = transformedTexture.texture;
								foundTextureCoordinates[0] = periodicity ? (bestTextureOrigin[0]+x)%bestTexture.width()
																		: bestTextureOrigin[0]+x;
								foundTextureCoordinates[1] = periodicity ? (bestTextureOrigin[1]+y)%bestTexture.height()
																		: bestTextureOrigin[1]+y;
								assert(transformedTexture.boundaries.pixelAbsolute(foundTextureCoordinates)==1.0); //TODO : shouldn't be needed
								if(	transformedTexture.boundaries.pixelAbsolute(foundTextureCoordinates)==1.0 && //TODO : same
									foundTextureCoordinates[0]>=0 && foundTextureCoordinates[0]<bestTexture.width()
									&& foundTextureCoordinates[1]>=0 && foundTextureCoordinates[1]<bestTexture.height())
									foundExpTexture.pixelAbsolute(tileCoordinates) = bestTexture.pixelAbsolute(foundTextureCoordinates);
								else
									foundExpTexture.pixelAbsolute(tileCoordinates) = ImageType::zero();
							}
						}
					std::cout << "emulated patch origin - " << emulatedPatchOrigin << std::endl;
					std::cout << "patch origin - " << patchOrigin << std::endl;
					std::cout << "output origin - " << outputOrigin << std::endl;

					patch.alphaMap().for_all_pixels([&] (const ImageAlphad::PixelType &pix, int x, int y)
					{
						if(pix>0)
						{
							PixelPos tileCoordinates;
							tileCoordinates[0] = (patchOrigin[0]+x)%emulatedAlpha.width();
							tileCoordinates[1] = (patchOrigin[1]+y)%emulatedAlpha.height();
							PixelPos outputCoordinates;
							outputCoordinates[0] = outputOrigin[0]+x;
							outputCoordinates[1] = outputOrigin[1]+y;
							ImageAlphad::PixelType alphaPix = emulatedAlpha.pixelAbsolute(tileCoordinates);
							if(alphaPix>0 && outputCoordinates[0]>0 && outputCoordinates[1]>0 && outputCoordinates[0]<output.width() && outputCoordinates[1]<output.height()) //<< used as a failsafe (sometimes a patch can take the entire image in size and this causes problems, might want to TODO into that...)
							{
								PixelType &pixOut = output.pixelAbsolute(outputCoordinates);
								DataType *pixOutData = reinterpret_cast<DataType *>(&pixOut);

								const PixelType &pixExpTexture = foundExpTexture.pixelAbsolute(tileCoordinates);
								const DataType *pixExpTextureData = reinterpret_cast<const DataType *>(&pixExpTexture);
								for(unsigned i=0; i<pixelSize; ++i)
								{
									uvMap.pixelAbsolute(outputCoordinates)[0]=double(foundTextureCoordinates[0])/foundExpTexture.width();
									uvMap.pixelAbsolute(outputCoordinates)[1]=double(foundTextureCoordinates[1])/foundExpTexture.height();
									if(rdosp != nullptr)
									{
										const DataType *targetNormData = reinterpret_cast<const DataType *>(&rdosp->targetNorm);
										const DataType *sourceNormData = reinterpret_cast<const DataType *>(&rdosp->sourceNorm);
										pixOutData[i] += pixExpTextureData[i] * alphaPix * (1.0/targetNormData[i]) * sourceNormData[i];
									}
									else
									{
										output.pixelAbsolute(outputCoordinates)
												+= foundExpTexture.pixelAbsolute(tileCoordinates) * alphaPix;
									}
								}
							}
						}
					});
				}
			}
	}
	// IO::save01_in_u8(uvMap, ${ASTEX_TEMPO_PATH}+"_uvMap" + std::to_string(time(nullptr)) + ".png");
	Mipmap<ImageType> mipmapOutput(output);
	mipmapOutput.setMaxPowReductionLevel(this->m_patchMapMipmap.maxPowReductionLevel());
	mipmapOutput.setMode(this->m_patchMapMipmap.mode());
	mipmapOutput.generate();
	return mipmapOutput;
}

}//namespace

}//namespace


#endif
