#ifndef __CTEXCH_CONTENT_H__
#define __CTEXCH_CONTENT_H__

#include <Eigen/Eigen>
#include "ASTex/mipmap.h"
#include "patch.h"


namespace ASTex
{

namespace ContentExchange
{

template<typename I>
class Patch;

class MipmapCEPatch;

template<typename I>
/**
 * @brief The MipmapContentExchange class is a particular type of mipmap which embodies mipmapped portions of an image.
 * These portions must have been stored in an image and mipmapped before using this class.
 * Content must have been shifted and emplaced where the patch's default content is in order to fit the size.
 * One should mostly use the above Content class to access it.
 */
class MipmapCEContent : public Mipmap<I>
{
public:
	/**
	 * @brief MipmapCEContent constructor for MipmapCEContent.
	 * @param contentColor input-sized mipmap containing the content, shifted to the position of the patch.
	 * @param patchAlpha input-sized mipmap containing the alpha information of the patch.
	 */
	MipmapCEContent();
	MipmapCEContent(const I& content);

	/**
	 * @brief generate generates the content mipmap with settings given by
	 * setTexture(), setMode() and setMaxPowReductionLevel().
	 * @pre setParentPatch must have been called, because patch informations are used to compute it.
	 * It is also expected to be a valid patch, that is, this content was build to fit the parentPatch.
	 */
	void generate();

	/**
	 * @brief setParentPatch provides this class with the parent patch.
	 * @param parentPatch
	 */
	void setParentPatch(const Patch<I> *parentPatch) {m_parentPatch=parentPatch;}

private:

	const Patch<I> *m_parentPatch;
};


template<typename I>
MipmapCEContent<I>::MipmapCEContent() :
	Mipmap<I>(),
	m_parentPatch(0)
{}

template<typename I>
MipmapCEContent<I>::MipmapCEContent(const I& content):
	Mipmap<I>(content),
	m_parentPatch(0)
{}

template<typename I>
void MipmapCEContent<I>::generate()
{
	static typename I::PixelType zero;
	assert(m_parentPatch
		   && "MipmapCEContent::generate: parent patch has not been given (try MipmapCEContent::setParentPatch)");
	assert(m_parentPatch->alphaMipmap().isGenerated()
		   && "MipmapCEContent::generate: parent patch's mipmap was not generated");
	assert(this->isTextureSet()
		   && "MipmapCEContent::generate: no content texture was set (try MipmapCEContent::setTexture)");

	//clean content before generating
	const MipmapCEPatch &patchMipmapAlpha = m_parentPatch->alphaMipmap();
	const ImageGrayd& correspondingPatchAlphaTexture = patchMipmapAlpha.texture();
	PixelPos patchOrigin = patchMipmapAlpha.originAt(0, 0);
	I cleanedTexture;
	cleanedTexture.initItk(this->m_isoMipmaps[0].width(), this->m_isoMipmaps[0].height(), true); //TODO: does not set at 0 (current bug)
	cleanedTexture.for_all_pixels([&] (typename I::PixelType &pix)
	{
		pix=zero;
	});
	correspondingPatchAlphaTexture.for_all_pixels([&] (const ImageGrayd::PixelType &pix, int x, int y)
	{
		if(pix>0)
		{
			int xShift=(x+patchOrigin[0])%cleanedTexture.width();
			int yShift=(y+patchOrigin[1])%cleanedTexture.height();
			if(!std::is_floating_point<typename I::DataType>::value)
			{
				size_t arraySize = sizeof(typename I::PixelType) / sizeof(typename I::DataType);
				typename I::DataType *pix2 = new typename I::DataType[arraySize];
				//uint64_t *pixi = new uint64_t[arraySize]();
				std::memcpy(pix2, &(this->m_isoMipmaps[0].pixelAbsolute(xShift, yShift)), sizeof(typename I::PixelType));
				for(size_t i=0; i<arraySize; ++i)
				{
					pix2[i] *= pix;
				}
				std::memcpy(&(cleanedTexture.pixelAbsolute(xShift, yShift)), pix2, sizeof(typename I::PixelType));
				delete[] pix2;
			}
			else
				cleanedTexture.pixelAbsolute(xShift, yShift) = this->m_isoMipmaps[0].pixelAbsolute(xShift, yShift) * pix;
		}
	});
	this->m_isoMipmaps[0].copy_pixels(cleanedTexture);

	Mipmap<I>::generate();

	unsigned i, j, maxIterations;
	unsigned xMin, yMin;
	auto computeMipmap = [&] ()
	{
		//find the size of the new, sparse content and its origin in contentColor's mipmap
		const ImageGrayd& correspondingPatchAlphaMipmap=patchMipmapAlpha.mipmap(i, j);
		I& oldContentColor = this->mipmap(i, j);
		I mipmapContentColor;
		mipmapContentColor.initItk(correspondingPatchAlphaMipmap.width(), correspondingPatchAlphaMipmap.height(), true);
		MipmapCEPatch::PixelPos origin = patchMipmapAlpha.originAt(i, j);
		xMin=(unsigned)origin[0];
		yMin=(unsigned)origin[1];

		//Fill up the new mipmap
		mipmapContentColor.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
		{
			pix = oldContentColor.pixelAbsolute(  (xMin+x)%oldContentColor.width(),
												  (yMin+y)%oldContentColor.height());
		});
		//replace the old with the new
		oldContentColor = mipmapContentColor;
	};

	//allocation

	maxIterations = std::max(this->numberMipmapsWidth(), this->numberMipmapsHeight());
	this->m_isoMipmaps.resize(maxIterations);
	if(this->mode() == ANISOTROPIC)
	{
		this->m_anisoMipmapsWidth.resize(this->numberMipmapsHeight()-1);
		for(i=0; i<this->m_anisoMipmapsWidth.size(); ++i)
		{
			this->m_anisoMipmapsWidth[i].resize(this->numberMipmapsWidth()-1-i); //TODO: check if correct when using rectangular textures
		}
		this->m_anisoMipmapsHeight.resize(this->numberMipmapsWidth()-1);
		for(i=0; i<this->m_anisoMipmapsHeight.size(); ++i)
		{
			this->m_anisoMipmapsHeight[i].resize(this->numberMipmapsHeight()-1-i); //TODO: same
		}
	}

	if(this->mode()==ISOTROPIC || this->mode()==NO_FILTER)
	for(i=0, j=0; i<maxIterations; ++i, ++j)
		computeMipmap();

	else if(this->mode()==ANISOTROPIC)
		for(i=0; i<this->numberMipmapsWidth(); ++i)
			for(j=0; j<this->numberMipmapsHeight(); ++j)
				computeMipmap();
}

template<typename I>
class Content
{
public:
	Content(const I &imageContainingContent, const Patch<I> &parentPatch);

	const I& texture() const;
	const I& mipmap(unsigned i, unsigned j) const;
	const MipmapCEContent<I>& contentMipmap() const {return m_explicitContentMipmap;}
	void setTranslationTag(Eigen::Vector2i translation) {m_translationTag = translation;}
	Eigen::Vector2i translationTag() const {return m_translationTag;}
private:
	MipmapCEContent<I> m_explicitContentMipmap;
	Eigen::Vector2i m_translationTag;
};

template<typename I>
Content<I>::Content(const I& imageContainingContent, const Patch<I> &parentPatch):
	m_explicitContentMipmap(imageContainingContent)
{
	m_explicitContentMipmap.setParentPatch(&parentPatch);
	m_explicitContentMipmap.setMode(parentPatch.alphaMipmap().mode());
	m_explicitContentMipmap.setMaxPowReductionLevel(parentPatch.alphaMipmap().maxPowReductionLevel());
	m_explicitContentMipmap.generate();
}

template<typename I>
const I& Content<I>::texture() const
{
	return m_explicitContentMipmap.texture();
}

template<typename I>
const I& Content<I>::mipmap(unsigned i, unsigned j) const
{
	return m_explicitContentMipmap.mipmap(i, j);
}

}

}

#endif
