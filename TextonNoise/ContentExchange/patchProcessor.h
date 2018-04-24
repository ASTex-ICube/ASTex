#ifndef __PATCH_PROCESSOR_H__
#define __PATCH_PROCESSOR_H__

#include "histogram.h"
#include "patch.h"
#include "mipmap.h"
#include "content.h"

namespace ASTex
{

namespace ContentExchange
{

using ImageMask64 = ImageGrayd;

template<typename I>
class PatchProcessor
{
public:
    PatchProcessor();
    PatchProcessor(const I& image);

    //get

    const I& getTexture() const;
    const Patch<I>& patchAt(size_t index) const;
    Patch<I> &patchAt(size_t index);

    size_t nbPatches() const {return m_patches.size();}

    //set

    void setTexture(const I& texture);
    void setFiltering(mipmap_mode_t mipmapMode);

    //misc

    void generate();

    void debug_setPatchFromImageRGBd(const ImageRGBd& patchImage);
    void debug_setRandomContents(unsigned nbContentsPerPatch);

    //iterators

    typedef typename std::vector<Patch<I>>::iterator iterator;
    typedef typename std::vector<Patch<I>>::const_iterator const_iterator;

    iterator begin() {return m_patches.begin();}
    const_iterator begin() const {return m_patches.begin();}

    iterator end() {return m_patches.end();}
    const_iterator end() const {return m_patches.end();}

//    iterator& operator++(iterator& it);
//    const_iterator& operator++(const_iterator& it) const;

//    iterator operator++(iterator& it, int);
//    const_iterator operator++(iterator& it, int) const;

    //static

    static void setDefaultFilteringMode(mipmap_mode_t defaultMipmapMode) {ms_defaultMipmapMode=defaultMipmapMode;}

private:

    std::vector<Patch<I>> m_patches;
    I m_texture;
    ImageMask64 m_patchMask; //allows for up to 64 patches. TODO: replace ImageMask64 with MipmapMask
    mipmap_mode_t m_mipmapMode;

    static mipmap_mode_t ms_defaultMipmapMode;
};

template<typename I>
mipmap_mode_t PatchProcessor<I>::ms_defaultMipmapMode = NO_FILTER;

template<typename I>
PatchProcessor<I>::PatchProcessor():
    m_patches(),
    m_texture(),
    m_patchMask(),
    m_mipmapMode(ms_defaultMipmapMode)
{}

template<typename I>
PatchProcessor<I>::PatchProcessor(const I& texture):
    m_patches(),
    m_texture(texture),
    m_patchMask(),
    m_mipmapMode(ms_defaultMipmapMode)
{}

//get

template<typename I>
const I& PatchProcessor<I>::getTexture() const
{
    return m_texture;
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

//set

template<typename I>
void PatchProcessor<I>::setTexture(const I& texture)
{
    m_texture=texture;
}

template<typename I>
void PatchProcessor<I>::setFiltering(mipmap_mode_t mipmapMode)
{
    m_mipmapMode = mipmapMode;
}

template<typename I>
void PatchProcessor<I>::generate()
{
    assert(m_texture.is_initialized() &&
           "PatchProcessor::generate: texture uninitialized (use PatchProcessor::setTexture with an initialized texture)");
    assert(m_patchMask.is_initialized() &&
           "PatchProcessor::generate: patch mask uninitialized (use PatchProcessor::computePatches to compute patches)");
    assert(m_texture.size()==m_patchMask.size() &&
           "PatchProcessor::generate: patch mask must have the same size as texture (texture changed?)");

    //find the number of patches

    //compute the bitwise maximum value of the texture, which looks like 0b00..011..11
    using word64=uint64_t;
    word64 w=0x0;
    m_patchMask.for_all_pixels([&] (ImageMask64::PixelType &pix)
    {
        w|=reinterpret_cast<word64&>(pix);
    });
    //its log's floor defines the number of different patches and is proprely computed this way:
    int lg;
    word64 wTest;
    for(lg=0, wTest=0x1; (w & wTest) > 0x0; ++lg, wTest*=2);
    m_patches.resize(lg);

    wTest=0x1;
    w=0x0;

    //then for each patch, compute their alpha map and give it to them. The patches will take care of the rest with reduce().
    for(typename std::vector<Patch<I>>::iterator it=m_patches.begin(); it!=m_patches.end(); ++it)
    {
        ImageGrayd alphaMap;
        alphaMap.initItk(m_texture.width(), m_texture.height(), true );
        m_patchMask.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
        {
            if((w|=reinterpret_cast<word64&>(pix))==wTest)
                alphaMap.pixelAbsolute(x, y)=1.0;
            w=0x0;
        });
        wTest*=2;
        Patch<I> &patch=(*it);
        patch.setAlphaMap(alphaMap, ms_defaultMipmapMode);
        Content<I> content(m_texture, (*it));
        patch.addContent(content);
    }
}

//misc

template<typename I>
void PatchProcessor<I>::debug_setPatchFromImageRGBd(const ImageRGBd& patchImage)
{
    assert(patchImage.size() == m_texture.size() &&
           "PatchProcessor::debug_setPatchFromImageRGBd: patchImage should have the same size as texture");
    HistogramRGBd histogramPI(patchImage);
    assert(histogramPI.binsNumber()<65 &&
           "PatchProcessor::debug_setPatchFromImageRGBd: there should be up to 64 different colors in the given image");
    using word=uint64_t;
    word w=0x1;
    m_patchMask.initItk(patchImage.width(), patchImage.height(), true);
    for(auto it=histogramPI.begin(); it!=histogramPI.end(); ++it)
    {
        patchImage.for_all_pixels([&] (const ImageRGBd::PixelType &pix, int x, int y)
        {
            if((*it).first==pix)
            {
                word &wPixel=reinterpret_cast<word&>(m_patchMask.pixelAbsolute(x, y));
                wPixel |= w;
            }
        });
        w*=2;
    }
}

template<typename I>
void PatchProcessor<I>::debug_setRandomContents(unsigned nbContentsPerPatch)
{

}


}

}


#endif
