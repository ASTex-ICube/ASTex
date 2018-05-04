#ifndef __CONTENT_H__
#define __CONTENT_H__

#include "mipmap.h"
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
 * @brief The MipmapContentExchange class is a particular type of mipmap
 * which embodies mipmapped portions of an image.
 * These portions must have been stored in an image and mipmapped before using this class.
 * Content must have been shifted and emplaced where the patch's default content is
 * in order to fit the size.
 */
class MipmapCEContent : public Mipmap<I>
{
public:
    /**
     * @brief MipmapCEContent constructor for MipmapCEContent.
     * Immediately produces the output image.
     * @param contentColor input-sized mipmap containing the content, shifted to the position of the patch.
     * @param patchAlpha input-sized mipmap containing the alpha information of the patch.
     */
    MipmapCEContent();
    MipmapCEContent(const I& content);

    void generate();

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
    if(!m_parentPatch)
    {
        std::cerr << "Warning: MipmapCEContent::generate: parent patch has not been given (try MipmapCEContent::setParentPatch)" << std::endl;
        return;
    }
    if(!this->isTextureSet())
    {
        std::cerr << "Warning: MipmapCEContent::generate: no content texture was set" << std::endl;
        return;
    }
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

    if(this->mode()==ISOTROPIC)
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
private:
    MipmapCEContent<I> m_explicitContentMipmap;
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
