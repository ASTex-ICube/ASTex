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

using ImageMask64 = ImageGrayu64;
using word64=uint64_t;

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

    void initializePatchesRegularGrid(unsigned nbPatches=64);
    void refinePatchesGI();
    void initializeContents();

    Mipmap<I> generate(int textureWidth, int textureHeight) const;

    void debug_setPatchFromImageRGBd(const ImageRGBd& patchImage);
    void debug_setRandomContents(unsigned nbContentsPerPatch);

    size_t analysis_getGPUMemoryCost() const;
    size_t analysis_getNumberOfTextureAccessForMipmap(unsigned i, unsigned j) const;

    //iterators

    typedef typename std::vector<Patch<I>>::iterator iterator;
    typedef typename std::vector<Patch<I>>::const_iterator const_iterator;

    iterator begin() {return m_patches.begin();}
    const_iterator begin() const {return m_patches.begin();}

    iterator end() {return m_patches.end();}
    const_iterator end() const {return m_patches.end();}

    //static

    static void setDefaultFilteringMode(mipmap_mode_t defaultMipmapMode) {ms_defaultMipmapMode=defaultMipmapMode;}

private:

    std::vector<Patch<I>> m_patches;
    I m_texture;
    MipmapBitmask<ImageMask64> m_patchMaskMipmap;
    mipmap_mode_t m_mipmapMode;

    static mipmap_mode_t ms_defaultMipmapMode;
};

template<typename I>
mipmap_mode_t PatchProcessor<I>::ms_defaultMipmapMode = NO_FILTER;

template<typename I>
PatchProcessor<I>::PatchProcessor():
    m_patches(),
    m_texture(),
    m_patchMaskMipmap(),
    m_mipmapMode(ms_defaultMipmapMode)
{}

template<typename I>
PatchProcessor<I>::PatchProcessor(const I& texture):
    m_patches(),
    m_texture(texture),
    m_patchMaskMipmap(),
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
void PatchProcessor<I>::initializePatchesRegularGrid(unsigned nbPatches)
{
    assert(m_texture.is_initialized() &&
           "PatchProcessor::initializePatches: texture uninitialized (use PatchProcessor::setTexture with an initialized texture)");
    assert(nbPatches <= 64 &&
           "PatchProcessor::initializePatches: synthesis cannot allow more than 64 patches (nor should need to)");
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

    m_patchMaskMipmap.setTexture(patchMap);
    m_patchMaskMipmap.setMode(ms_defaultMipmapMode);
    m_patchMaskMipmap.generate();
}

template<typename I>
void PatchProcessor<I>::refinePatchesGI()
{
    assert(m_texture.is_initialized() &&
           "PatchProcessor::refinePatchesGI: texture uninitialized (use PatchProcessor::setTexture with an initialized texture)");
    assert(m_patchMaskMipmap.isGenerated() &&
           "PatchProcessor::refinePatchesGI: patch mask mipmap not generated (use PatchProcessor::initializePatches<Mode> to compute patches)");
    assert(m_texture.size()==m_patchMaskMipmap.texture().size() &&
           "PatchProcessor::refinePatchesGI: patch mask must have the same size as texture (texture changed?)");



}

template<typename I>
void PatchProcessor<I>::initializeContents()
{
    assert(m_texture.is_initialized() &&
           "PatchProcessor::initializeContents: texture uninitialized (use PatchProcessor::setTexture with an initialized texture)");
    assert(m_patchMaskMipmap.isGenerated() &&
           "PatchProcessor::initializeContents: patch mask mipmap not generated (use PatchProcessor::initializePatches<Mode> to compute patches)");
    assert(m_texture.size()==m_patchMaskMipmap.texture().size() &&
           "PatchProcessor::initializeContents: patch mask must have the same size as texture (texture changed?)");

    //find the number of patches

    //compute the bitwise maximum value of the texture, which looks like 0b00..011..11
    word64 w=0x0;
    m_patchMaskMipmap.texture().for_all_pixels([&] (ImageMask64::PixelType &pix)
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
        m_patchMaskMipmap.texture().for_all_pixels([&] (ImageMask64::PixelType &pix, int x, int y)
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

template<typename I>
Mipmap<I> PatchProcessor<I>::generate(int textureWidth, int textureHeight) const
{
    assert(m_patches.size()>0 && "PatchProcessor::generate: initialize() must be called after setting patch mask");
    I output;
    output.initItk(m_texture.width(), m_texture.height(), true);
    Mipmap<I> mipmapOutput(output);
    mipmapOutput.setMaxPowReductionLevel(m_patchMaskMipmap.maxPowReductionLevel());
    mipmapOutput.setMode(m_patchMaskMipmap.mode());
    mipmapOutput.generate();
    for(unsigned i=0; i<mipmapOutput.numberMipmapsWidth(); ++i)
        for(unsigned j=0; j<mipmapOutput.numberMipmapsHeight(); ++j)
        {
            I &mipmap = mipmapOutput.mipmap(i, j);
            mipmap.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
            {
                pix=typename I::PixelType();
                //the following utilizes bitmasks to read the data of each mipmap
                word64 wPixel=reinterpret_cast<word64>(m_patchMaskMipmap.mipmap(i, j).pixelAbsolute(x, y));
                word64 w=0x1;
                for(size_t p=0; p<m_patches.size(); ++p)
                {
                    word64 wTmp=w;
                    if((wTmp&=wPixel)>0)
                    {
                        PixelPos patchOrigin = this->patchAt(p).alphaMipmap().originAt(i, j);
                        int x2=x-patchOrigin[0];
                        int y2=y-patchOrigin[1];
                        if(x2 < 0)
                            x2+=mipmap.width();
                        if(y2 < 0)
                            y2+=mipmap.height();

                        //choose your content here
                        pix += this->patchAt(p).contentAt(0).mipmap(i, j).pixelAbsolute(x2, y2);
                    }
                    w*=2;
                }
            });
        }

    return mipmapOutput;
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
    word64 w=0x1;
    ImageMask64 patchMask;
    patchMask.initItk(patchImage.width(), patchImage.height(), true);
    for(auto it=histogramPI.begin(); it!=histogramPI.end(); ++it)
    {
        patchImage.for_all_pixels([&] (const ImageRGBd::PixelType &pix, int x, int y)
        {
            if((*it).first==pix)
            {
                word64 &wPixel=reinterpret_cast<word64&>(patchMask.pixelAbsolute(x, y));
                wPixel |= w;
            }
        });
        w*=2;
    }
    m_patchMaskMipmap.setTexture(patchMask);
    m_patchMaskMipmap.setMode(ms_defaultMipmapMode);
    m_patchMaskMipmap.generate();
}

template<typename I>
void PatchProcessor<I>::debug_setRandomContents(unsigned nbContentsPerPatch)
{
    assert(m_patches.size()>0 && "PatchProcessor::debug_setRandomContents: initialize() must be called before being able to chose contents");
    unsigned i, j;
    int randomShiftX, randomShiftY;
    I shiftedTexture;
    shiftedTexture.initItk(m_texture.width(), m_texture.height());

    for(i=0; i<nbContentsPerPatch; ++i)
    {
        randomShiftX = rand();
        randomShiftY = rand();
        shiftedTexture.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
        {
            pix = m_texture.pixelAbsolute((x + randomShiftX)%shiftedTexture.width(), (y + randomShiftY)%shiftedTexture.height());
        });
        for(j=0; j<m_patches.size(); ++j)
        {
            Patch<I> &patch=m_patches[j];
            Content<I> c(shiftedTexture, patch);
            patch.addContent(c);
        }
    }

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
    for(i=0; i<m_patchMaskMipmap.numberMipmapsWidth(); ++i)
    {
        if(m_patchMaskMipmap.mode()==ISOTROPIC) //if isotropic, there are mipmaps only on the diagonal (i, i)
            addImageMaskMemoryCostOf(m_patchMaskMipmap.mipmap(i, i));
        else
            for(j=0; j<m_patchMaskMipmap.numberMipmapsHeight(); ++j)
                addImageMaskMemoryCostOf(m_patchMaskMipmap.mipmap(i, j));
    }
    for(i=0; i<m_patches.size(); ++i)
    {
        const Patch<I> &patch=m_patches[i];
        for(j=0; j<patch.numberContents(); ++j)
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
size_t PatchProcessor<I>::analysis_getNumberOfTextureAccessForMipmap(unsigned i, unsigned j) const
{
    const ImageMask64& mipmap=m_patchMaskMipmap.mipmap(i, j);
    unsigned access = 0; //counts the number of texture access
    mipmap.for_all_pixels([&] (const ImageMask64::PixelType &pix)
    {
        //mask analysis
        word64 wPixel=reinterpret_cast<word64>(pix);
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

}

}


#endif
