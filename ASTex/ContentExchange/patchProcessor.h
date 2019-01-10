#ifndef __CTEXCH_PATCH_PROCESSOR_H__
#define __CTEXCH_PATCH_PROCESSOR_H__

#include "ASTex/histogram.h"
#include "patch.h"
#include "content.h"
#include "ASTex/easy_io.h"
#include "ASTex/Stamping/sampler.h"
#include "old_PatchProcessor.h"
#include "old_PoissonDiskSampler.h"
#include <random>
#include "ASTex/PCTS/pcts.h"

namespace ASTex
{

namespace ContentExchange
{

#define CTEXCH_MAX_PATCHES 64
using ImageMask64 = ImageGrayu64;
using word64=uint64_t;

template<typename I>
class Atlas;

template<typename I>
class AlternativeAtlas;

template<typename I>
/**
 * @brief The PatchProcessor class is a class made for processing and storing content exchange elements.
 * Feed it a texture, choose how you want patches and contents to be computed and bake an offline result with
 * generate(), or export the class for online rendering.
 */
class PatchProcessor
{
public:
	PatchProcessor();
	PatchProcessor(const I& image);

	~PatchProcessor();

	//get

	const I& texture() const;

	/**
	 * @brief filteringMode
	 * @return the filtering method, between ISOTROPIC, ANISOTROPIC and NO_FILTER.
	 */
	mipmap_mode_t filteringMode() const {return m_mipmapMode;}

	/**
	 * @brief patchAt
	 * @param index
	 * @return the patch at index index.
	 */
	const Patch<I>& patchAt(size_t index) const;
	Patch<I> &patchAt(size_t index);

	const I& atlasContentAt(size_t index) const;
	I &atlasContentAt(size_t index);

	const Atlas<I> &atlas() const {return *m_atlas;}

	const Eigen::Vector2i contentTranslationAt(unsigned p, unsigned c) const;

	/**
	 * @brief patchMapMipmap corresponds to the structure used for storing a patch map for each lod.
	 * Each patchmap is an image in 64 bits where each bit corresponds to a patch ID.
	 * You can use shortcut functions such as mipmap() or patchMap() to get directly to sub-structures of the mipmap.
	 * @return the patchmap's mipmap.
	 */
	const MipmapBitmask<ImageMask64>& patchMapMipmap() const {return m_patchMapMipmap;}

	/**
	 * @brief mipmap
	 * @param k
	 * @param l
	 * @return patchmapMipmap().mipmap(k, l), which corresponds to the ImageMask64 mipmap of reduction level (k, l).
	 */
	const ImageMask64& mipmap(unsigned k, unsigned l) const {return m_patchMapMipmap.mipmap(k, l);}

	/**
	 * @brief patchmap
	 * @return patchmapMipmap().texture(), which corresponds to the highest level of ImageMask64,
	 * as well as the only image there is under the NO_FILTER filtering mode.
	 */
	const ImageMask64& patchmap() const {return m_patchMapMipmap.texture();}

	/**
	 * @brief maxMipmapPowReductionLevel
	 * Note that each of the 3 types of mipmaps
	 * (MipmapCEPatch of each patch, MipmapCEContent of each content, and this MipmapBitmask)
	 * are expected to have the same maximum reduction level, making this function alright to bound each of them.
	 * @return the maximum reduction level computed, in power of two.
	 */
	unsigned maxMipmapPowReductionLevel() const {return m_patchMapMipmap.maxPowReductionLevel(); }

	/**
	 * @brief numberMipmapsWidth
	 * Note that each of the 3 types of mipmaps
	 * (MipmapCEPatch of each patch, MipmapCEContent of each content, and this MipmapBitmask)
	 * are expected to have the same maximum reduction level, making this function alright to bound each of them.
	 * @return the number of mipmaps reduced in width. Returns 1 under NO_FILTER condition or previous to generate().
	 * Returns the same as numberMipmapsHeight() under ISOTROPIC conditions. Returns 0 if no texture was set.
	 */
	unsigned numberMipmapsWidth() const {return m_patchMapMipmap.numberMipmapsWidth(); }

	/**
	 * @brief numberMipmapsHeight
	 * Note that each of the 3 types of mipmaps
	 * (MipmapCEPatch of each patch, MipmapCEContent of each content, and this MipmapBitmask)
	 * are expected to have the same maximum reduction level, making this function alright to bound each of them.
	 * @return the number of mipmaps reduced in height. Returns 1 under NO_FILTER condition or previous to generate().
	 * Returns the same as numberMipmapsWidth() under ISOTROPIC conditions. Returns 0 if no texture was set.
	 */
	unsigned numberMipmapsHeight() const {return m_patchMapMipmap.numberMipmapsHeight(); }


	/**
	 * @brief nbPatches
	 * @return the number of patches (currently).
	 */
	size_t nbPatches() const {return m_patches.size();}

	/**
	 * @brief nbContents
	 * @return the expected number of contents. For boundary check, patches' content vector's sizes should be checked.
	 */
	size_t nbContents() const {return m_nbContentsPerPatch;}

	void setOutputSize(size_t width, size_t height) {m_outputWidth = width; m_outputHeight = height;}

	void setNbContentsPerPatch(unsigned nbContents) {m_nbContentsPerPatch = nbContents;}

	void setTexture(const I& texture);

	/**
	 * @brief setFilteringMode changes the filtering mode, between NO_FILTER, ISOTROPIC and ANISOTROPIC.
	 * @pre While not detected, DO NOT call this function with a different parameter between the generation
	 * (or exporting to GPU step) and the initialization, as this would result in undefined behavior.
	 * @param defaultMipmapMode
	 */
	void setFilteringMode(mipmap_mode_t defaultMipmapMode);

	void setSeed(unsigned seed);

	//////////////////////////////
	/////    INITIALIZERS    /////
	//////////////////////////////
	/**
	 * @brief initializePatchesRegularGrid builds patches from a regular grid.
	 * @param nbPatches is the expected number of patches. Will be cropped to the closest integer squared root.
	 * @pre nbPatches <= CTEXCH_MAX_PATCHES
	 */
	void patches_initRegularGrid(unsigned nbPatches=CTEXCH_MAX_PATCHES);

	/**
	 * @brief initializePatchesPoisson builds patches from a poisson sample.
	 * @param nbPatches is the number of patches.
	 * @pre nbPatches <= CTEXCH_MAX_PATCHES
	 */
	void patches_initCircles(unsigned nbPatches);

	void patches_initRandom(unsigned nbPatches);

	template<typename R>
	/**
							* @brief debug_setPatchFromImageRGB computes patches from an ImageRGB called patchImage.
							* each color in patchImage represents a patch.
							* @param patchImage
							*/
	void patches_initFromImageRGB(const ImageCommon<ImageRGBBase<R>, false> &patchImage);

	template<typename R>
	void patches_initFromImageGray(const ImageCommon<ImageGrayBase<R>, false> &patchImage);

	void patches_dilate(bool con8=false);

	void patches_fill(unsigned nbPatches);

	/**
	 * @brief initializeContents builds the initial contents from the input.
	 * Mandatory, but can be used in conjonction with other contents_ functions.
	 */
	void contents_initDefault();

	/**
	 * @brief debug_setRandomContents chooses and builds nbContentsPerPatch contents for each patches.
	 * Contents are taken from the input. They are chosen with a non-tweaked iid law.
	 * @param nbContentsPerPatch
	 */
	void contents_initRandom();

	void contents_enhancePCTS(std::string pctsArgFile = std::string());

	void fullProcess_GIOptimization(unsigned int nbPatchesPerDimension=5);

	template<typename std::enable_if<std::is_same<I, ASTex::ImageRGBu8>::value>::type* = nullptr >
	void fullProcess_oldMethod();

	//////////////////////////////
	/////   MISC FUNCTIONS   /////
	//////////////////////////////

	/**
	* @brief generate returns a synthesized output texture + mipmaps (same size as the input).
	* @param textureWidth expected output width. (TODO)
	* @param textureHeight expected output height.
	* @return the output mipmap. Use generate(...).texture() or generate(...).mipmap(0, 0) to get the output texture.
	*/
	Mipmap<I> generate() const;

	/**
	 * @brief saveRenderingPack saves most of the datas in outputDirectory.
	 * This includes every content of every patch for every level of detail computed, and other
	 * datas necesssary for the loading. This is a WIP and there is no load function. (TODO)
	 * @param outputDirectory
	 */
	void saveRenderingPack(const std::string &outputDirectory, bool storeLevel0=true) const;

	void loadRenderingPack(const std::string &inputDirectory);

	/**
	 * @brief analysis_getGPUMemoryCost
	 * TODO: could change easily.
	 * @return the overall memory cost of the datas to be stored in the GPU.
	 */
	size_t analysis_getGPUMemoryCost() const;

	/**
	 * @brief analysis_getNumberOfTextureAccessForMipmap
	 * @param i reduction in width
	 * @param j reduction in height
	 * @return the number of accesses needed to compute an output mipmap of reduction level k, l.
	 */
	size_t analysis_getNumberOfTextureAccessForMipmap(unsigned k, unsigned l) const;

	bool patchesOverlap() const {return m_patchesOverlap;}

	//iterators

	/**
	 * @brief iterator and const_iterator iterate over the patches.
	 */
	typedef typename std::vector<Patch<I>>::iterator iterator;
	typedef typename std::vector<Patch<I>>::const_iterator const_iterator;

	iterator begin() {return m_patches.begin();}
	const_iterator begin() const {return m_patches.begin();}

	iterator end() {return m_patches.end();}
	const_iterator end() const {return m_patches.end();}

private:

	unsigned m_nbContentsPerPatch;
	std::vector<Patch<I>> m_patches;
	I m_texture;
	MipmapBitmask<ImageMask64> m_patchMapMipmap;
	mipmap_mode_t m_mipmapMode;
	size_t m_outputWidth;
	size_t m_outputHeight;
	unsigned m_seed;

	//For rendering only
	Atlas<I> *m_atlas;
	std::vector<I> m_contentsAtlas;
	std::vector<Eigen::Vector2i> m_contentTranslations;

	bool m_patchesOverlap;

	static typename I::PixelType ms_zero;
};

template<typename I>
typename I::PixelType PatchProcessor<I>::ms_zero;

template<typename I>
PatchProcessor<I>::PatchProcessor():
	m_nbContentsPerPatch(1),
	m_patches(),
	m_texture(),
	m_patchMapMipmap(),
	m_mipmapMode(NO_FILTER),
	m_outputWidth(512),
	m_outputHeight(512),
	m_seed(0),
	m_atlas(0),
	m_contentsAtlas(),
	m_patchesOverlap(true)
{}

template<typename I>
PatchProcessor<I>::PatchProcessor(const I& texture):
	m_nbContentsPerPatch(1),
	m_patches(),
	m_texture(texture),
	m_patchMapMipmap(),
	m_mipmapMode(NO_FILTER),
	m_outputWidth(512),
	m_outputHeight(512),
	m_seed(0),
	m_atlas(0),
	m_contentsAtlas(),
	m_patchesOverlap(true)
{}

template<typename I>
PatchProcessor<I>::~PatchProcessor()
{
	delete m_atlas;
}

//get

template<typename I>
const I& PatchProcessor<I>::texture() const
{
	return m_texture;
}

template<typename I>
const Eigen::Vector2i PatchProcessor<I>::contentTranslationAt(unsigned p, unsigned c) const
{
	if(c == 0)
		return Eigen::Vector2i(0, 0);
	else
		return m_contentTranslations[p*m_nbContentsPerPatch+c];
}

template<typename I>
const Patch<I> &PatchProcessor<I>::patchAt(size_t index) const
{
	return m_patches[index];
}

template<typename I>
Patch<I>& PatchProcessor<I>::patchAt(size_t index)
{
	return m_patches[index];
}

template<typename I>
const I& PatchProcessor<I>::atlasContentAt(size_t index) const
{
	return m_contentsAtlas[index];
}

template<typename I>
I& PatchProcessor<I>::atlasContentAt(size_t index)
{
	return m_contentsAtlas[index];
}

//set

template<typename I>
void PatchProcessor<I>::setTexture(const I& texture)
{
	m_texture = texture;
}

template<typename I>
void PatchProcessor<I>::setFilteringMode(mipmap_mode_t mipmapMode)
{
	m_mipmapMode = mipmapMode;
}

template<typename I>
void PatchProcessor<I>::setSeed(unsigned seed)
{
	m_seed = seed;
}

//Patches initializers

template<typename I>
void PatchProcessor<I>::patches_initRegularGrid(unsigned nbPatches)
{
	assert(m_texture.is_initialized() &&
		   "PatchProcessor::initializePatches: texture uninitialized (use PatchProcessor::setTexture with an initialized texture)");
	assert(nbPatches <= CTEXCH_MAX_PATCHES &&
		   "PatchProcessor::initializePatches: synthesis cannot allow that many patches (max given by CTEXCH_MAX_PATCHES)");
	srand(m_seed);
	double nbPixelsInWidthPerPatch, nbPixelsInHeightPerPatch;
	int nbPatchesPerSide = int(std::sqrt(nbPatches));
	nbPatches = nbPatchesPerSide * nbPatchesPerSide;
	ImageMask64 patchMap;
	patchMap.initItk(m_texture.width(), m_texture.height());
	nbPixelsInWidthPerPatch = m_texture.width()/double(nbPatchesPerSide);
	nbPixelsInHeightPerPatch = m_texture.height()/double(nbPatchesPerSide);

	patchMap.for_all_pixels([&] (ImageMask64::PixelType &pix, int x, int y)
	{
		int xId = x/nbPixelsInWidthPerPatch;
		int yId = y/nbPixelsInHeightPerPatch;
		int id = yId * nbPatchesPerSide + xId;
		pix = uint64_t(std::pow(2, id));
	});

	m_patchMapMipmap.setTexture(patchMap);
	m_patchMapMipmap.setMode(m_mipmapMode);
	m_patchMapMipmap.generate();
}

template<typename I>
void PatchProcessor<I>::patches_initCircles(unsigned nbPatches)
{
	assert(m_texture.is_initialized() &&
		   "PatchProcessor::initializePatches: texture uninitialized (use PatchProcessor::setTexture with an initialized texture)");
	assert(nbPatches <= CTEXCH_MAX_PATCHES &&
		   "PatchProcessor::initializePatches: synthesis cannot allow that many patches (max given by CTEXCH_MAX_PATCHES)");
	srand(m_seed);
	Stamping::SamplerUniform sampler(nbPatches);
	std::vector<Eigen::Vector2f> centroids = sampler.generate();
	ImageMask64 patchMap;
	patchMap.initItk(m_texture.width(), m_texture.height());
	float r = m_texture.width() * m_texture.height() / float(nbPatches*nbPatches*2);
	int i = 0;
	for(std::vector<Eigen::Vector2f>::const_iterator cit = centroids.begin(); cit != centroids.end(); ++cit, ++i)
	{
		word64 w=0x1;
		w <<= i;
		patchMap.for_all_pixels([&] (ImageMask64::PixelType &pix, int x, int y)
		{
			int x0, y0;
			x0 = (*cit)[0]*m_texture.width();
			y0 = (*cit)[1]*m_texture.height();
			if(std::sqrt((x0-x) * (x0-x) + (y0-y) * (y0-y)) < r)
				reinterpret_cast<word64&>(pix) |= w;
		});
	}

	m_patchMapMipmap.setTexture(patchMap);
	m_patchMapMipmap.setMode(m_mipmapMode);
	m_patchMapMipmap.generate();
}

template<typename I>
void PatchProcessor<I>::patches_initRandom(unsigned nbPatches)
{
	//patch mask initialization
	assert(nbPatches < 64);
	static ImageMask64::PixelType zero;
	ImageMask64 patchMap;
	patchMap.initItk(m_texture.width(), m_texture.height(), true);
	patchMap.for_all_pixels([&] (ImageMask64::PixelType &pix)
	{
		pix = zero;
	});

	//sampler poisson but with the correct number of samples
	std::vector< Eigen::Vector2d > samples;
	std::vector< Eigen::Vector2d > samplesValidated;
	double minDist = std::sqrt(m_texture.width()*m_texture.height())/(2.0*std::sqrt(nbPatches) + 1);
	ContentExchg::PoissonDiskSampler::generateSamples(m_texture.width(), m_texture.height(), samples, minDist, nbPatches);
	std::minstd_rand minstdRandGenerator;
	minstdRandGenerator.seed(m_seed+1);
	assert(samples.size()>=nbPatches);
	for(size_t i=0; i<nbPatches; ++i)
	{
		std::uniform_int_distribution<size_t> distribution(0, samples.size()-1);
		size_t k = distribution(minstdRandGenerator);
		samplesValidated.push_back(samples[k]);
		samples.erase(samples.begin()+int(k));
	}

	//region growing
	unsigned *grownTexelCountdownVector = new unsigned[nbPatches]();
	unsigned totalNumberOfMarkedPixels = 0;
	for(unsigned i=0; i<nbPatches; ++i)
	{
		unsigned maxTexelsForPatchi = ((0.5+(rand()/double(RAND_MAX)))*m_texture.width()*m_texture.height())/nbPatches;
		grownTexelCountdownVector[i] = maxTexelsForPatchi;
	}
	do
	{
		word64 w = 0x1;
		int i = 0;
		for(auto &seed : samplesValidated)
		{
			ASTex::RegionGrowingOperator< ASTex::NeighborhoodPeriodic > regionOp( patchMap );
			ImageMask64::PixelType seedTex;
			seedTex = w;
			itk::Index<2> seedInt;
			seedInt[0]=int(seed[0]);
			seedInt[1]=int(seed[1]);
			regionOp.newRegionFromSeed< ASTex::GrowingStrategyRandom >(
						seedInt,

						// Condition function for accepting a pixel to the current region.
						[&]( const ASTex::NeighborhoodPeriodic::PixelType &texelId )
			{
				return (patchMap.pixelAbsolute(texelId) == 0x0 || patchMap.pixelAbsolute(texelId)==seedTex)
						&& grownTexelCountdownVector[i] > 0;
			},

			// Processing function applied to accepted pixels.
			[&]( ASTex::NeighborhoodPeriodic::PixelType &texelId, const ASTex::NeighborhoodPeriodic::NeighborsList& /*neighbors*/ )
			{
				ImageMask64::PixelType &pix = patchMap.pixelAbsolute(texelId);
				if(pix == 0x0)
				{
					--grownTexelCountdownVector[i];
					++totalNumberOfMarkedPixels;
					pix = seedTex;
				}
			});
			w*=2;
			grownTexelCountdownVector[i++]+=64; //whatever average number is fine, this is used to add more pixels in case the image isn't filled
		}
	}
	while(totalNumberOfMarkedPixels < m_texture.width()*m_texture.height());

	//Debug only
	ASTex::ImageGrayu8 patches;
	patches.initItk(m_texture.width(), m_texture.height());
	word64 w=0x1;
	for(unsigned i=0; i<nbPatches; ++i)
	{
		patchMap.for_all_pixels([&] (const word64 &pix, int x, int y)
		{
			if(pix == w)
				patches.pixelAbsolute(x, y)=uint8_t((i*255)/(nbPatches-1));
		});
		w *= 2;
	}
	patches_initFromImageGray(patches);

	delete[] grownTexelCountdownVector;
}

template<typename I>
template<typename R>
void PatchProcessor<I>::patches_initFromImageRGB(const ImageCommon<ImageRGBBase<R>, false>& patchImage)
{
	assert(patchImage.size() == m_texture.size() &&
		   "PatchProcessor::debug_setPatchFromImageRGBd: patchImage should have the same size as texture");
	HistogramRGBBase<R> histogramPI(patchImage);
	assert(histogramPI.binsNumber()<=CTEXCH_MAX_PATCHES &&
		   "PatchProcessor::debug_setPatchFromImageRGBd: there should be up to 64 different colors in the given image");
	word64 w=0x1;
	ImageMask64 patchMask;
	patchMask.initItk(patchImage.width(), patchImage.height(), true);
	for(auto it=histogramPI.begin(); it!=histogramPI.end(); ++it)
	{
		patchImage.for_all_pixels([&] (const typename ImageCommon<ImageRGBBase<R>, false>::PixelType &pix, int x, int y)
		{
			if((*it).first==pix)
			{
				word64 &wPixel=reinterpret_cast<word64&>(patchMask.pixelAbsolute(x, y));
				wPixel |= w;
			}
		});
		w*=2;
	}
	m_patchMapMipmap.setTexture(patchMask);
	m_patchMapMipmap.setMode(m_mipmapMode);
	m_patchMapMipmap.generate();
}

template<typename I>
template<typename R>
void PatchProcessor<I>::patches_initFromImageGray(const ImageCommon<ImageGrayBase<R>, false> &patchImage)
{
	assert(patchImage.size() == m_texture.size() &&
		   "PatchProcessor::debug_setPatchFromImageRGBd: patchImage should have the same size as texture");
	HistogramGrayBase<R> histogramPI(patchImage);
	assert(histogramPI.binsNumber()<=CTEXCH_MAX_PATCHES &&
		   "PatchProcessor::debug_setPatchFromImageRGBd: there should be up to 64 different colors in the given image");
	word64 w=0x1;
	ImageMask64 patchMask;
	patchMask.initItk(patchImage.width(), patchImage.height(), true);
	for(auto it=histogramPI.begin(); it!=histogramPI.end(); ++it)
	{
		patchImage.for_all_pixels([&] (const typename ImageCommon<ImageGrayBase<R>, false>::PixelType &pix, int x, int y)
		{
			if((*it).first==pix)
			{
				word64 &wPixel=reinterpret_cast<word64&>(patchMask.pixelAbsolute(x, y));
				wPixel |= w;
			}
		});
		w*=2;
	}
	m_patchMapMipmap.setTexture(patchMask);
	m_patchMapMipmap.setMode(m_mipmapMode);
	m_patchMapMipmap.generate();
}

template<typename I>
void PatchProcessor<I>::patches_dilate(bool con8)
{
	m_patchesOverlap = true;
	ImageMask64 originalPatchmap = m_patchMapMipmap.texture();
	int width = originalPatchmap.width();
	int height = originalPatchmap.height();
	ImageMask64 newPatchmap;
	newPatchmap.initItk(width, height);
	newPatchmap.copy_pixels(m_patchMapMipmap.texture());
	word64 w = 0x1;
	bool existsPatchw=true;

	for(int i=0; existsPatchw && i<64; ++i, w*=2)
	{
		existsPatchw = false;
		originalPatchmap.for_all_pixels([&](const ImageMask64::PixelType &pix, int x, int y)
		{
			if(pix == w)
			{
				x+=width;
				y+=height;
				existsPatchw = true;
//                newPatchmap.pixelAbsolute(x%width, (y-1)%height) |= w;
				newPatchmap.pixelAbsolute(x%width, (y+1)%height) |= w;
//                newPatchmap.pixelAbsolute((x-1)%width, y%height) |= w;
				newPatchmap.pixelAbsolute((x+1)%width, y%height) |= w;
				if(con8)
				{
					newPatchmap.pixelAbsolute((x-1)%width, (y-1)%height) |= w;
					newPatchmap.pixelAbsolute((x-1)%width, (y+1)%height) |= w;
					newPatchmap.pixelAbsolute((x+1)%width, (y-1)%height) |= w;
					newPatchmap.pixelAbsolute((x+1)%width, (y+1)%height) |= w;
				}
			}
		});
	}
	m_patchMapMipmap.setTexture(newPatchmap);
	m_patchMapMipmap.setMode(m_mipmapMode);
	m_patchMapMipmap.generate();
}

template<typename I>
void PatchProcessor<I>::patches_fill(unsigned nbPatches)
{
	m_patchesOverlap = true;
	static ImageMask64::PixelType zero;
	ImageMask64 patchMap;
	word64 w=0x1, wFill=0x0;
	for(unsigned i=0; i<nbPatches; ++i)
	{
		wFill |= w;
		w *= 2;
	}
	patchMap.initItk(m_texture.width(), m_texture.height(), true);
	patchMap.for_all_pixels([&] (ImageMask64::PixelType &pix)
	{
		pix = wFill;
	});
	m_patchMapMipmap.setTexture(patchMap);
	m_patchMapMipmap.setMode(m_mipmapMode);
	m_patchMapMipmap.generate();
}

//Content initializers

template<typename I>
void PatchProcessor<I>::contents_initDefault()
{
	assert(m_texture.is_initialized() &&
		   "PatchProcessor::initializeContents: texture uninitialized (use PatchProcessor::setTexture with an initialized texture)");
	assert(m_patchMapMipmap.isGenerated() &&
		   "PatchProcessor::initializeContents: patch mask mipmap not generated (use PatchProcessor::initializePatches<Mode> to compute patches)");
	assert(m_texture.size()==m_patchMapMipmap.texture().size() &&
		   "PatchProcessor::initializeContents: patch mask must have the same size as texture (texture changed?)");

	//find the number of patches

	//compute the bitwise maximum value of the texture, which looks like 0b00..011..11
	word64 w=0x0;
	m_patchMapMipmap.texture().for_all_pixels([&] (ImageMask64::PixelType &pix)
	{
		w|=reinterpret_cast<word64&>(pix);
	});
	//its log's floor defines the number of different patches and is proprely computed this way:
	int lg;
	word64 wTest;
	//we can then predict the number of patches based on what we read in the texture
	for(lg=0, wTest=0x1; (w & wTest) > 0x0; ++lg, wTest*=2);
	m_patches.resize(lg);

	wTest=0x1;
	w=0x0;

	//then for each patch, compute their alpha map and give it to them. The patches will take care of the rest.
	for(typename std::vector<Patch<I>>::iterator it=m_patches.begin(); it!=m_patches.end(); ++it)
	{
		ImageGrayd alphaMap;
		alphaMap.initItk(m_texture.width(), m_texture.height(), true );
		m_patchMapMipmap.texture().for_all_pixels([&] (ImageMask64::PixelType &pix, int x, int y)
		{
			if((w|=reinterpret_cast<word64&>(pix)) & wTest)
			{
				//check how many patches there are
				int n=0;
				for(int i=0, w = 0x1; i<lg; ++i, w*=2)
				{
					if((w & pix)!=0)
						++n;
				}
				alphaMap.pixelAbsolute(x, y)=1.0/n;
			}
			w=0x0;
		});
		wTest*=2;
		Patch<I> &patch=(*it);
		patch.setAlphaMap(alphaMap, m_mipmapMode);
		Content<I> content(m_texture, (*it));
		patch.addContent(content);
	}
}

template<typename I>
void PatchProcessor<I>::contents_initRandom()
{
	assert(m_patches.size()>0 &&
		   "PatchProcessor::contents_initRandom: a patch initializer must be called before being able to choose contents");
	srand(m_seed);
	unsigned i, j;
	int randomShiftX, randomShiftY;
	I shiftedTexture;
	shiftedTexture.initItk(m_texture.width(), m_texture.height());
	for(j=0; j<m_patches.size(); ++j)
	{
		Patch<I> &patch=m_patches[j];
		for(i=0; i<m_nbContentsPerPatch-1; ++i)
		{
			randomShiftX = rand()%shiftedTexture.width();
			randomShiftY = rand()%shiftedTexture.height();
			shiftedTexture.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
			{
				pix = m_texture.pixelAbsolute((x + randomShiftX)%shiftedTexture.width(), (y + randomShiftY)%shiftedTexture.height());
			});
			Content<I> c(shiftedTexture, patch);
			c.setTranslationTag(Eigen::Vector2i(randomShiftX, randomShiftY));
			patch.addContent(c);
		}
	}
}

template<typename I>
void PatchProcessor<I>::contents_enhancePCTS(std::string pctsArgFile)
{
	assert(m_patches.size()>0 && m_patches[0].nbContents()>1 &&
		   "PatchProcessor::debug_setRandomContents: alternative contents must be present to optimize them.");
	size_t arraySize = sizeof(typename I::PixelType) / sizeof(typename I::DataType);
	typename I::DataType *pix_one = new typename I::DataType[arraySize];
	static typename I::PixelType pix_zero;
	ASTex::Pcts<I> pcts(m_texture);
	pcts.setWidth(m_texture.width());
	pcts.setHeight(m_texture.height());
	pcts.setPeriodicity(true);
	if(!pctsArgFile.empty())
		pcts.loadParametersFromFile(pctsArgFile);
	for(int i=0; i<m_patches.size(); ++i)
	{
		Patch<I> &p = m_patches[i];
		MipmapCEPatch::PixelPos origin = p.originAt(0, 0);
		I mask;
		mask.initItk(m_texture.width(), m_texture.height());
		mask.for_all_pixels([] (typename I::PixelType &pix)
		{
			pix = pix_zero;
		});
		for(int x=0; x<p.alphaMap().width(); ++x)
		{
			for(int y=0; y<p.alphaMap().height(); ++y)
			{
				if(p.alphaMap().pixelAbsolute(x, y)>0)
				{
					mask.pixelAbsolute( (x+origin[0])%mask.width(),
										(y+origin[1])%mask.height())
							= pix_one;
				}
			}
		}
//        mask.save(std::string("/home/nlutz/ieee2019/_debug/mask_p") + std::to_string(i) + ".png");
		for(int j=1; j<p.nbContents(); ++j)
		{
			Content<I> &c = p.contentAt(j);
			I contentTexture;
			contentTexture.initItk(m_texture.width(), m_texture.height());
			contentTexture.copy_pixels(m_texture);

			for(int x=0; x<c.texture().width(); ++x)
			{
				for(int y=0; y<c.texture().height(); ++y)
				{
					if(p.alphaMap().pixelAbsolute(x, y)>0)
					{
						contentTexture.pixelAbsolute(   (x+origin[0])%contentTexture.width(),
														(y+origin[1])%contentTexture.height())
								= c.texture().pixelAbsolute(x, y);
					}
				}
			}
//            contentTexture.save(std::string("/home/nlutz/ieee2019/_debug/inputWithContent_p") + std::to_string(i) + "_c" + std::to_string(j) + ".png");
			pcts.setSynthesis(contentTexture, mask);
			contentTexture = pcts.generate();
//            contentTexture.save(std::string("/home/nlutz/ieee2019/_debug/pctsd_p") + std::to_string(i) + "_c" + std::to_string(j) + ".png");
			Content<I> pctsdContent(contentTexture, p);
			c = pctsdContent;
		}
	}
	delete[](pix_one);
}

//Full process initializers (patches + contents)

template<typename I>
template<typename std::enable_if<std::is_same<I, ASTex::ImageRGBu8>::value>::type*>
void PatchProcessor<I>::fullProcess_oldMethod()
{
	unsigned fragmentMinSize        = 20;
	unsigned fragmentMaxSize        = 500;
	unsigned fragmentColorThreshold = 40;
	unsigned requiredPatchNumber    = 16;
	unsigned downsamplingMinSize    = 128;
	srand(m_seed);

	ContentExchg::FragmentProcessor fProc( m_texture );
	fProc.createFragments( fragmentMaxSize, fragmentColorThreshold );
	fProc.cleanupSmallestFragments( fragmentMinSize );
	fProc.updateFragmentsAttributes();

	// Patches generation from old content exchange

	ContentExchg::PatchProcessor pProc( fProc );
	pProc.createPatches( int(requiredPatchNumber) );
	pProc.computePatchBoundaries();

	ASTex::ImageRGBu8::PixelType *idColors = new ASTex::ImageRGBu8::PixelType [ fProc.fragmentCount() ];
	for( int i=0; i<pProc.patchCount(); ++i )
	{
		idColors[i].SetRed  ( uint8_t((i * 255)/pProc.patchCount()-1) );
		idColors[i].SetGreen( uint8_t((i * 255)/pProc.patchCount()-1) );
		idColors[i].SetBlue ( uint8_t((i * 255)/pProc.patchCount()-1) );
	}

	ASTex::ImageRGBu8 *imgPatch = new ASTex::ImageRGBu8( m_texture.width(), m_texture.height() );
	for( auto &p : pProc.patches() )
		for( auto &frag : p.fragments )
			for( auto &pixel : fProc.fragmentById(int(frag)).pixels )
				imgPatch->pixelAbsolute(pixel) = idColors[p.id];

	//content generation from old content exchange

	std::vector<double> rotations;
	rotations.push_back( 0.0 );

	std::vector<double> scales;
	scales.push_back( 1.0000 );

	pProc.findAlternativeContents( rotations, scales, m_nbContentsPerPatch-1, downsamplingMinSize );

	//translation to new content exchange

	patches_initFromImageRGB(*imgPatch);
	patches_dilate();
	contents_initDefault();
	int pId = 0;
	for(auto &p : pProc.patches())
	{
		for(auto &c : p.contents)
		{
			I shiftedTexture;
			shiftedTexture.initItk(m_texture.width(), m_texture.height());

			Patch<I> &patch=m_patches[pId];
			shiftedTexture.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
			{
				pix = m_texture.pixelAbsolute(  (x+c.offset[0])%shiftedTexture.width(),
						(y+c.offset[1])%shiftedTexture.height());
			});
			Content<I> content(shiftedTexture, patch);
			patch.addContent(content);
			shiftedTexture.initItk(m_texture.width(), m_texture.height());
		}
		++pId;
	}

	//Debug only
	for( int i=0; i<pProc.patchCount(); ++i )
	{
		idColors[i].SetRed  ( uint8_t((i * 255)/pProc.patchCount()-1) );
		idColors[i].SetGreen( uint8_t((i * 255)/pProc.patchCount()-1) );
		idColors[i].SetBlue ( uint8_t((i * 255)/pProc.patchCount()-1) );
	}
	for( auto &p : pProc.patches() )
		for( auto &frag : p.fragments )
			for( auto &pixel : fProc.fragmentById(int(frag)).pixels )
				imgPatch->pixelAbsolute(pixel) = idColors[p.id];
}

template<typename I>
void PatchProcessor<I>::fullProcess_GIOptimization(unsigned int nbPatchesPerDimension)
{
	ImageMask64 patchMask;
	patchMask.initItk(m_texture.width(), m_texture.height(), true);
	patchMask.for_all_pixels([&] (const typename I::PixelType &pix)
	{
		word64 &wPixel=reinterpret_cast<word64&>(pix);
		wPixel = 0x0;
	});

}

//Generator

template<typename I>
Mipmap<I> PatchProcessor<I>::generate() const
{
	assert(m_patches.size()>0 && "PatchProcessor::generate: a patch initializer must be called before using generate()");
	I output;
	output.initItk(m_outputWidth, m_outputHeight, true);
	Mipmap<I> mipmapOutput(output);
	mipmapOutput.setMaxPowReductionLevel(m_patchMapMipmap.maxPowReductionLevel());
	mipmapOutput.setMode(m_patchMapMipmap.mode());
	mipmapOutput.generate();

	srand(m_seed);
	int generationSeed = rand();

	auto lmbd_generate1Mipmap = [&] (int k, int l)
	{
		I &mipmap = mipmapOutput.mipmap(k, l);
		mipmap.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
		{
			pix=PatchProcessor<I>::ms_zero;
			//the following utilizes bitmasks to read the data of each mipmap
			int xClamped = x%m_patchMapMipmap.mipmap(k, l).width();
			int yClamped = y%m_patchMapMipmap.mipmap(k, l).height();
			word64 wPixel=word64( m_patchMapMipmap.mipmap(k, l).pixelAbsolute(
									 xClamped,
									 yClamped) );
			word64 w=0x1;
			for(size_t p=0; p<m_patches.size(); ++p)
			{
				word64 wTmp=w;
				if((wTmp&=wPixel)>0)
				{
					PixelPos patchOrigin = this->patchAt(p).alphaMipmap().originAt(k, l);
					int x2=xClamped-patchOrigin[0];
					int y2=yClamped-patchOrigin[1];
					if(x2 < 0)
						x2+=m_patchMapMipmap.mipmap(k, l).width();
					if(y2 < 0)
						y2+=m_patchMapMipmap.mipmap(k, l).height();

					//choose your content here
					srand(generationSeed%65535 + 65535*std::floor((x - patchOrigin[0])/m_patchMapMipmap.mipmap(k, l).width())
						  + 1  *std::floor((y - patchOrigin[1])/m_patchMapMipmap.mipmap(k, l).height()));
					pix += this->patchAt(p).contentAt(rand()%this->patchAt(p).nbContents()).mipmap(k, l).pixelAbsolute(x2, y2);
				}
				w*=2;
			}
		});
	};

	for(unsigned i=0; i<mipmapOutput.numberMipmapsWidth(); ++i)
	{
		if(mipmapOutput.mode()==ANISOTROPIC)
		{
			for(unsigned j=0; j<mipmapOutput.numberMipmapsHeight(); ++j)
			{
				lmbd_generate1Mipmap(i, j);
			}
		}
		else
		{
			lmbd_generate1Mipmap(i, i);
		}
	}

	return mipmapOutput;
}

//misc

template<typename I>
void PatchProcessor<I>::saveRenderingPack(const std::string &outputDirectory, bool storeLevel0) const
{
	unsigned i,j,k,l;
	//saving input
	m_texture.save(outputDirectory + "/input.png");

	//saving patch map
	for(k=0; k<m_patchMapMipmap.numberMipmapsWidth(); ++k)
	{
		if(m_patchMapMipmap.mode()==ANISOTROPIC)
			for(l=0; l<m_patchMapMipmap.numberMipmapsHeight(); ++l)
			{
				ImageMask64 im64 = m_patchMapMipmap.mipmap(k, l);
				Histogram<ImageMask64>::saveImageToCsv(im64, outputDirectory + "/patchmap_mw" + std::to_string(k) + "_mh" + std::to_string(l) + ".csv");
			}
		else
		{
			ImageMask64 im64 = m_patchMapMipmap.mipmap(k, k);
			Histogram<ImageMask64>::saveImageToCsv(im64, outputDirectory + "/patchmap_mw" + std::to_string(k) + "_mh" + std::to_string(k) + ".csv");
		}
	}

	Atlas<I> atlas(*this);
	atlas.setStoreHighestLevel(storeLevel0);
	atlas.generateAndSaveAtlas(outputDirectory);

	//saving useful data
	std::ofstream ofs_data_out(outputDirectory + "/data.csv");
	ofs_data_out << m_patches.size() << std::endl;
	ofs_data_out << nbContents() << std::endl;
	ofs_data_out << (unsigned int)m_patchMapMipmap.mode() << std::endl;
	ofs_data_out << m_patchMapMipmap.numberMipmapsWidth() << std::endl;
	ofs_data_out << m_patchMapMipmap.numberMipmapsHeight() << std::endl;
	ofs_data_out << m_patchMapMipmap.maxPowReductionLevel() << std::endl;
	ofs_data_out.close();

	//saving origins
	for(i=0; i<m_patches.size(); ++i)
	{
		std::ofstream ofs_origins_out(outputDirectory + "/origin_p" + std::to_string(i) + ".csv");
		const Patch<I> &patch = m_patches[i];

		for(k=0; k<patch.alphaMipmap().numberMipmapsWidth(); ++k)
		{
			if(patch.alphaMipmap().mode()==ANISOTROPIC)
				for(l=0; l<patch.alphaMipmap().numberMipmapsHeight(); ++l)
				{
					ofs_origins_out << patch.alphaMipmap().originAt(k, l)[0]
									<< " "
									<< patch.alphaMipmap().originAt(k, l)[1];
				}
			else
			{
				ofs_origins_out << patch.alphaMipmap().originAt(k, k)[0]
								<< " "
								<< patch.alphaMipmap().originAt(k, k)[1];
			}
			ofs_origins_out << std::endl;
		}
		ofs_origins_out.close();
	}
	std::ofstream ofs_executable(outputDirectory + "/pack.exch");
	ofs_executable << outputDirectory << std::endl;
	ofs_executable.close();

	//additional: content translations
	std::ofstream ofs_transformations_out(outputDirectory + "/transformations.csv");
	for(i = 0; i<m_patches.size(); ++i)
	{
		const Patch<I> &patch = m_patches[i];
		for(j = 1; j<nbContents(); ++j)
		{
			ofs_transformations_out << patch.contentAt(j).translationTag()[0] << ' '
									<< patch.contentAt(j).translationTag()[1] << std::endl;
		}
	}
	ofs_transformations_out.close();
}

template<typename I>
void PatchProcessor<I>::loadRenderingPack(const std::string &inputDirectory)
{
	unsigned i, j, k,l;

	//loading input
	m_texture.load(inputDirectory + "/input.png");

	//loading useful data
	std::ifstream ifs_data_in(inputDirectory + "/data.csv");

	size_t nbPatches, nbContents, numberMipmapsWidth, numberMipmapsHeight, maxPowReductionLevel, ui_mipmapMode;
	mipmap_mode_t m_mipmapMode;

	ifs_data_in >> nbPatches;
	ifs_data_in >> nbContents;
	ifs_data_in >> ui_mipmapMode;
	ifs_data_in >> numberMipmapsWidth;
	ifs_data_in >> numberMipmapsHeight;
	ifs_data_in >> maxPowReductionLevel;
	ifs_data_in.close();

	m_mipmapMode = mipmap_mode_t(ui_mipmapMode);

	m_patches.resize(nbPatches);
	setNbContentsPerPatch(nbContents);
	for(i=0; i<m_patches.size(); ++i)
	{
		std::ifstream ifs_origins_in(inputDirectory + "/origin_p" + std::to_string(i) + ".csv");
		Patch<I> &patch = m_patches[i];
		patch.alphaMipmap().setOriginsSize(numberMipmapsWidth, numberMipmapsHeight);

		for(k=0; k<numberMipmapsWidth; ++k)
		{
			int s, t;
			if(m_mipmapMode==ANISOTROPIC)
				for(l=0; l<numberMipmapsHeight; ++l)
				{
					ifs_origins_in >> s >> t;
					patch.alphaMipmap().originAt(k, l)[0] = s;
					patch.alphaMipmap().originAt(k, l)[1] = t;
				}
			else
			{
				ifs_origins_in >> s >> t;
				patch.alphaMipmap().originAt(k, k)[0] = s;
				patch.alphaMipmap().originAt(k, k)[1] = t;
			}
		}
		ifs_origins_in.close();
	}

	//loading patch map
	ImageMask64 patchMap;
	Histogram<ImageMask64>::loadImageFromCsv(patchMap, inputDirectory + "/patchmap_mw" + std::to_string(0)
																				+ "_mh" + std::to_string(0) + ".csv");
	m_patchMapMipmap.setTexture(patchMap);
	m_patchMapMipmap.setMode(m_mipmapMode);
	m_patchMapMipmap.setMaxPowReductionLevel(maxPowReductionLevel);
	m_patchMapMipmap.generate();
	for(k=0; k<numberMipmapsWidth; ++k)
	{
		if(m_mipmapMode==ANISOTROPIC)
			for(l=0; l<numberMipmapsHeight; ++l)
			{
				Histogram<ImageMask64>::loadImageFromCsv(m_patchMapMipmap.mipmap(k, l), inputDirectory + "/patchmap_mw" + std::to_string(k) + "_mh" + std::to_string(l) + ".csv");
			}
		else
		{
			Histogram<ImageMask64>::loadImageFromCsv(m_patchMapMipmap.mipmap(k, k), inputDirectory + "/patchmap_mw" + std::to_string(k) + "_mh" + std::to_string(k) + ".csv");
		}
	}

	m_atlas = new Atlas<I>(*this);
	m_atlas->loadOrigins(inputDirectory + "/atlasOrigins.csv");
	m_contentsAtlas.resize(nbContents);
	for(i=0; i<nbContents; ++i)
	{
		m_contentsAtlas[i].load(inputDirectory + "/contentAtlas_" + std::to_string(i) + ".png");
	}

	ifs_data_in.open(inputDirectory + "/transformations.csv");
	m_contentTranslations.resize(0);
	for(i=0; i<nbPatches; ++i)
	{
		m_contentTranslations.push_back(Eigen::Vector2i(0, 0));
		for(j=1; j<nbContents; ++j)
		{
			int x, y;
			ifs_data_in >> x >> y;
			m_contentTranslations.push_back(Eigen::Vector2i(x, y));
		}
	}
	ifs_data_in.close();

}

template<typename I>
size_t PatchProcessor<I>::analysis_getGPUMemoryCost() const
{
	unsigned i, j, k, l;
	size_t s=0;

	//lambda to compute the memory cost of an image
	auto addImageMemoryCostOf = [&](const I &image)
	{
		s+=sizeof(typename I::PixelType)*image.height()*image.width();
	};
	//lambda to compute the memory cost of the mask image
	auto addImageMaskMemoryCostOf = [&](const ImageMask64 &image)
	{
		s+=sizeof(typename ImageMask64::PixelType)*image.height()*image.width();
	};
	for(i=0; i<m_patchMapMipmap.numberMipmapsWidth(); ++i)
	{
		if(m_patchMapMipmap.mode()==ISOTROPIC) //if isotropic, there are mipmaps only on the diagonal (i, i)
			addImageMaskMemoryCostOf(m_patchMapMipmap.mipmap(i, i));
		else
			for(j=0; j<m_patchMapMipmap.numberMipmapsHeight(); ++j)
				addImageMaskMemoryCostOf(m_patchMapMipmap.mipmap(i, j));
	}
	for(i=0; i<m_patches.size(); ++i)
	{
		const Patch<I> &patch=m_patches[i];
		for(j=0; j<patch.nbContents(); ++j)
		{
			const Content<I> &content=patch.contentAt(j);
			for(k=0; k<content.contentMipmap().numberMipmapsWidth(); ++k)
			{
				if(content.contentMipmap().mode()==ISOTROPIC)
					addImageMemoryCostOf(content.mipmap(k, k));
				else
					for(l=0; l<content.contentMipmap().numberMipmapsHeight(); ++l)
						addImageMemoryCostOf(content.mipmap(k, l));
			}
		}
	}
	return s;
}

template<typename I>
size_t PatchProcessor<I>::analysis_getNumberOfTextureAccessForMipmap(unsigned k, unsigned l) const
{
	const ImageMask64& mipmap=m_patchMapMipmap.mipmap(k, l);
	unsigned access = 0; //counts the number of texture access
	mipmap.for_all_pixels([&] (const ImageMask64::PixelType &pix)
	{
		//mask analysis
		word64 wPixel=word64(pix);
		word64 w=0x1;
		for(size_t p=0; p<m_patches.size(); ++p)
		{
			word64 wTmp=w;
			if((wTmp&=wPixel)>0)
				++access;
			w *= 2;
		}
	});

	return access;
}

}//namespace

}//namespace


#endif
