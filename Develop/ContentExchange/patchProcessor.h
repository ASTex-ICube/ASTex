#ifndef __CTEXCH_PATCH_PROCESSOR_H__
#define __CTEXCH_PATCH_PROCESSOR_H__

#include "histogram.h"
#include "patch.h"
#include "mipmap.h"
#include "content.h"
#include "ASTex/easy_io.h"
#include "Stamping/sampler.h"

namespace ASTex
{

namespace ContentExchange
{

#define CTEXCH_MAX_PATCHES 64
using ImageMask64 = ImageGrayu64;
using word64=uint64_t;

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

    /**
     * @brief patchMapMipmap corresponds to the structure used for storing a patch map for each lod.
     * Each patchmap is an image in 64 bits where each bit corresponds to a patch ID.
     * You can use shortcut functions such as mipmap() or patchMap() to get directly to sub-structures of the mipmap.
     * @return the patchmap's mipmap.
     */
    const MipmapBitmask<ImageMask64>& patchMapMipmap() const {return m_patchMaskMipmap;}

    /**
     * @brief mipmap
     * @param k
     * @param l
     * @return patchmapMipmap().mipmap(k, l), which corresponds to the ImageMask64 mipmap of reduction level (k, l).
     */
    const ImageMask64& mipmap(unsigned k, unsigned l) const {return m_patchMaskMipmap.mipmap(k, l);}

    /**
     * @brief patchmap
     * @return patchmapMipmap().texture(), which corresponds to the highest level of ImageMask64,
     * as well as the only image there is under the NO_FILTER filtering mode.
     */
    const ImageMask64& patchmap() const {return m_patchMaskMipmap.texture();}

    /**
     * @brief maxMipmapPowReductionLevel
     * Note that each of the 3 types of mipmaps
     * (MipmapCEPatch of each patch, MipmapCEContent of each content, and this MipmapBitmask)
     * are expected to have the same maximum reduction level, making this function alright to bound each of them.
     * @return the maximum reduction level computed, in power of two.
     */
    unsigned maxMipmapPowReductionLevel() const {return m_patchMaskMipmap.maxPowReductionLevel(); }

    /**
     * @brief numberMipmapsWidth
     * Note that each of the 3 types of mipmaps
     * (MipmapCEPatch of each patch, MipmapCEContent of each content, and this MipmapBitmask)
     * are expected to have the same maximum reduction level, making this function alright to bound each of them.
     * @return the number of mipmaps reduced in width. Returns 1 under NO_FILTER condition or previous to generate().
     * Returns the same as numberMipmapsHeight() under ISOTROPIC conditions. Returns 0 if no texture was set.
     */
    unsigned numberMipmapsWidth() const {return m_patchMaskMipmap.numberMipmapsWidth(); }

    /**
     * @brief numberMipmapsHeight
     * Note that each of the 3 types of mipmaps
     * (MipmapCEPatch of each patch, MipmapCEContent of each content, and this MipmapBitmask)
     * are expected to have the same maximum reduction level, making this function alright to bound each of them.
     * @return the number of mipmaps reduced in height. Returns 1 under NO_FILTER condition or previous to generate().
     * Returns the same as numberMipmapsWidth() under ISOTROPIC conditions. Returns 0 if no texture was set.
     */
    unsigned numberMipmapsHeight() const {return m_patchMaskMipmap.numberMipmapsHeight(); }


    /**
     * @brief nbPatches
     * @return the number of patches (currently).
     */
    size_t nbPatches() const {return m_patches.size();}

    size_t nbContents() const {return m_patches.size()>0 ? m_patches[0].nbContents() : 0;}

    //set

    /**
     * @brief setTexture sets the input texture.
     * @param texture
     */
    void setTexture(const I& texture);

    /**
     * @brief setFilteringMode changes the filtering mode, between NO_FILTER, ISOTROPIC and ANISOTROPIC.
     * @pre While not detected, DO NOT call this function with a different parameter between the generation
     * (or exporting to GPU step) and the initialization, as this would result in undefined behavior.
     * @param defaultMipmapMode
     */
    void setFilteringMode(mipmap_mode_t defaultMipmapMode);

                            //////////////////////////////
                            /////    INITIALIZERS    /////
                            //////////////////////////////
    /**
     * @brief initializePatchesRegularGrid builds patches from a regular grid.
     * @param nbPatches is the expected number of patches. Will be cropped to the closest integer squared root.
     * @pre nbPatches <= CTEXCH_MAX_PATCHES
     */
    void initializePatchesRegularGrid(unsigned nbPatches=CTEXCH_MAX_PATCHES);

    /**
     * @brief initializePatchesPoisson builds patches from a poisson sample.
     * @param nbPatches is the number of patches.
     * @pre nbPatches <= CTEXCH_MAX_PATCHES
     */
    void initializePatchesPoissonCircles(unsigned nbPatches);

    /**
     * @brief initializeContents builds the initial contents from the input.
     */
    void initializeContents();

    template<typename R>
    /**
     * @brief debug_setPatchFromImageRGB computes patches from an ImageRGB called patchImage.
     * each color in patchImage represents a patch.
     * @param patchImage
     */
    void debug_setPatchFromImageRGB(const ImageCommon<ImageRGBBase<R>, false> &patchImage);

    /**
     * @brief debug_setRandomContents chooses and builds nbContentsPerPatch contents for each patches.
     * Contents are taken from the input. They are chosen with a non-tweaked iid law.
     * @param nbContentsPerPatch
     */
    void debug_setRandomContents(unsigned nbContentsPerPatch, unsigned int seed=0);


                           //////////////////////////////
                           /////   MISC FUNCTIONS   /////
                           //////////////////////////////

    /**
    * @brief generate returns a synthesized output texture + mipmaps (same size as the input).
    * @param textureWidth expected output width. (TODO)
    * @param textureHeight expected output height.
    * @return the output mipmap. Use generate(...).texture() or generate(...).mipmap(0, 0) to get the output texture.
    */
    Mipmap<I> generate(int textureWidth, int textureHeight) const;


    /**
     * @brief saveRenderingPack saves most of the datas in outputDirectory.
     * This includes every content of every patch for every level of detail computed, and other
     * datas necesssary for the loading. This is a WIP and there is no load function. (TODO)
     * @param outputDirectory
     */
    void saveRenderingPack(const std::string &outputDirectory);

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

    std::vector<Patch<I>> m_patches;
    I m_texture;
    MipmapBitmask<ImageMask64> m_patchMaskMipmap;
    mipmap_mode_t m_mipmapMode;

    static typename I::PixelType ms_zero;
};

template<typename I>
typename I::PixelType PatchProcessor<I>::ms_zero;

template<typename I>
PatchProcessor<I>::PatchProcessor():
    m_patches(),
    m_texture(),
    m_patchMaskMipmap(),
    m_mipmapMode(NO_FILTER)
{}

template<typename I>
PatchProcessor<I>::PatchProcessor(const I& texture):
    m_patches(),
    m_texture(texture),
    m_patchMaskMipmap(),
    m_mipmapMode(NO_FILTER)
{}

//get

template<typename I>
const I& PatchProcessor<I>::texture() const
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
void PatchProcessor<I>::setFilteringMode(mipmap_mode_t mipmapMode)
{
    m_mipmapMode = mipmapMode;
}

template<typename I>
void PatchProcessor<I>::initializePatchesRegularGrid(unsigned nbPatches)
{
    assert(m_texture.is_initialized() &&
           "PatchProcessor::initializePatches: texture uninitialized (use PatchProcessor::setTexture with an initialized texture)");
    assert(nbPatches <= CTEXCH_MAX_PATCHES &&
           "PatchProcessor::initializePatches: synthesis cannot allow that many patches (max given by CTEXCH_MAX_PATCHES)");
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
    m_patchMaskMipmap.setMode(m_mipmapMode);
    m_patchMaskMipmap.generate();
}

template<typename I>
void PatchProcessor<I>::initializePatchesPoissonCircles(unsigned nbPatches)
{
    assert(m_texture.is_initialized() &&
           "PatchProcessor::initializePatches: texture uninitialized (use PatchProcessor::setTexture with an initialized texture)");
    assert(nbPatches <= CTEXCH_MAX_PATCHES &&
           "PatchProcessor::initializePatches: synthesis cannot allow that many patches (max given by CTEXCH_MAX_PATCHES)");
    Stamping::SamplerPoisson sampler(nbPatches);
    std::vector<Eigen::Vector2f> centroids = sampler.generate();
    ImageMask64 patchMap;
    patchMap.initItk(m_texture.width(), m_texture.height());
    float r = m_texture.width() * m_texture.height() / float(nbPatches*nbPatches);
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
            if((w|=reinterpret_cast<word64&>(pix)) & wTest)
                alphaMap.pixelAbsolute(x, y)=1.0;
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
Mipmap<I> PatchProcessor<I>::generate(int textureWidth, int textureHeight) const
{
    assert(m_patches.size()>0 && "PatchProcessor::generate: initialize() must be called after setting patch mask");
    I output;
    output.initItk(m_texture.width(), m_texture.height(), true);
    Mipmap<I> mipmapOutput(output);
    mipmapOutput.setMaxPowReductionLevel(m_patchMaskMipmap.maxPowReductionLevel());
    mipmapOutput.setMode(m_patchMaskMipmap.mode());
    mipmapOutput.generate();

    auto lmbd_generate1Mipmap = [&] (int k, int l)
    {
        I &mipmap = mipmapOutput.mipmap(k, l);
        mipmap.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
        {
            pix=PatchProcessor<I>::ms_zero;
            //the following utilizes bitmasks to read the data of each mipmap
            word64 wPixel=reinterpret_cast<word64>(m_patchMaskMipmap.mipmap(k, l).pixelAbsolute(x, y));
            word64 w=0x1;
            for(size_t p=0; p<m_patches.size(); ++p)
            {
                word64 wTmp=w;
                if((wTmp&=wPixel)>0)
                {
                    PixelPos patchOrigin = this->patchAt(p).alphaMipmap().originAt(k, l);
                    int x2=x-patchOrigin[0];
                    int y2=y-patchOrigin[1];
                    if(x2 < 0)
                        x2+=mipmap.width();
                    if(y2 < 0)
                        y2+=mipmap.height();

                    //choose your content here
                    pix += this->patchAt(p).contentAt(0).mipmap(k, l).pixelAbsolute(x2, y2);
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
template <typename R>
void PatchProcessor<I>::debug_setPatchFromImageRGB(const ImageCommon<ImageRGBBase<R>, false>& patchImage)
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

    m_patchMaskMipmap.setTexture(patchMask);
    m_patchMaskMipmap.setMode(m_mipmapMode);
    m_patchMaskMipmap.generate();
}

template<typename I>
void PatchProcessor<I>::debug_setRandomContents(unsigned nbContentsPerPatch, unsigned int seed)
{
    assert(m_patches.size()>0 && "PatchProcessor::debug_setRandomContents: initialize() must be called before being able to chose contents");
    srand(seed);
    unsigned i, j;
    int randomShiftX, randomShiftY;
    I shiftedTexture;
    shiftedTexture.initItk(m_texture.width(), m_texture.height());
    for(j=0; j<m_patches.size(); ++j)
    {
        Patch<I> &patch=m_patches[j];
        for(i=0; i<nbContentsPerPatch; ++i)
        {
            randomShiftX = rand();
            randomShiftY = rand();
            shiftedTexture.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
            {
                pix = m_texture.pixelAbsolute((x + randomShiftX)%shiftedTexture.width(), (y + randomShiftY)%shiftedTexture.height());
            });
            Content<I> c(shiftedTexture, patch);
            patch.addContent(c);
        }
    }

}

template<typename I>
void PatchProcessor<I>::saveRenderingPack(const std::string &outputDirectory)
{
    unsigned i,j,k,l;
    //saving patch map
    for(k=0; k<m_patchMaskMipmap.numberMipmapsWidth(); ++k)
    {
        if(m_patchMaskMipmap.mode()==ANISOTROPIC)
            for(l=0; l<m_patchMaskMipmap.numberMipmapsHeight(); ++l)
            {
                ImageMask64 im64 = m_patchMaskMipmap.mipmap(k, l);
                Histogram<ImageMask64>::saveImageToCsv(im64, outputDirectory + "/patchmap_mw" + std::to_string(k) + "_mh" + std::to_string(l) + ".csv");
            }
        else
        {
            ImageMask64 im64 = m_patchMaskMipmap.mipmap(k, k);
            Histogram<ImageMask64>::saveImageToCsv(im64, outputDirectory + "/patchmap_mw" + std::to_string(k) + "_mh" + std::to_string(k) + ".csv");
        }

    }
    //saving contents
    for(i=0; i<m_patches.size(); ++i)
    {
        const Patch<I> &patch = m_patches[i];
        for(j=0; j<patch.nbContents(); ++j)
        {
            const Content<I> content = patch.contentAt(j);
            for(k=0; k<content.contentMipmap().numberMipmapsWidth(); ++k)
            {
                if(content.contentMipmap().mode() == ANISOTROPIC)
                {
                    for(l=0; l<content.contentMipmap().numberMipmapsHeight(); ++l)
                    {
                        I im = content.mipmap(k, l);
                        IO::save01_in_u8(im, outputDirectory + "/p" + std::to_string(i) + "_c" + std::to_string(j) + "_mw" + std::to_string(k) + "_mh" + std::to_string(l) + ".png");
                    }
                }
                else
                {
                    I im = content.mipmap(k, k);
                    IO::save01_in_u8(im, outputDirectory + "/p" + std::to_string(i) + "_c" + std::to_string(j) + "_mw" + std::to_string(k) + "_mh" + std::to_string(k) + ".png");
                }
            }
        }
    }
    //saving useful data
    std::ofstream ofs_data_out(outputDirectory + "/data.csv");
    ofs_data_out << m_patches.size() << std::endl;
    ofs_data_out << m_patches[0].nbContents() << std::endl;
    ofs_data_out << m_patchMaskMipmap.mode() << std::endl;
    ofs_data_out << m_patchMaskMipmap.numberMipmapsWidth() << std::endl;
    ofs_data_out << m_patchMaskMipmap.numberMipmapsHeight() << std::endl;
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
                    ofs_origins_out << patch.alphaMipmap().originAt(k, l) << " ";
                }
            else
            {
                ofs_origins_out << patch.alphaMipmap().originAt(k, k);
            }
            ofs_origins_out << std::endl;

        }

        ofs_origins_out.close();
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
    const ImageMask64& mipmap=m_patchMaskMipmap.mipmap(k, l);
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

}//namespace

}//namespace


#endif
