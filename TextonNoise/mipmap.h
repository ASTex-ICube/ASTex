#ifndef __MIPMAP__H__
#define __MIPMAP__H__

#include <cmath>
#include <ASTex/special_io.h>
#include <ASTex/utils.h>
#include <ASTex/colorspace_filters.h>
#include <ASTex/mask.h>
#include <type_traits>

using namespace ASTex;

typedef enum {ISOTROPIC=0, ANISOTROPIC=1, NO_FILTER=2} mipmap_mode_t;

template <class I>
/**
 * @brief The Mipmap class is a base class for mip-mapping.
 * It is able to filter any image with mipmaps by computing a recursive bilinear average over a 4 pixels neighboorhood.
 * the filterDivide functions determine how the mip-map is handled, and they are the functions to override to change
 * the behaviour of the mipmapping.
 * The default mode for Mipmaps is ISOTROPIC.
 */
class Mipmap
{
public:

    Mipmap();
    Mipmap(const I& texture);

    //get

    mipmap_mode_t mode() const {return m_mode;}
    size_t numberMipmapsWidth() const;
    size_t numberMipmapsHeight() const;
    const I& mipmap(unsigned xPowReduction, unsigned yPowReduction) const;
    const I& texture() const;
    I& texture();

    bool generated() const {return m_generated;}
    bool textureSet() const {return m_textureSet;}

    //set

    void setTexture(const I& texture);

    //misc

    virtual void generate(mipmap_mode_t mode=ISOTROPIC, unsigned maxPowReductionLevel=0);
    void revertSetTexture();
    void revertGenerate();

    void fullMipmap(I& fullMipmap);

protected:

    virtual void filterDivide2Width(const I& texture, I& result);
    virtual void filterDivide2Height(const I& texture, I& result);
    virtual void filterDivide2Full(const I& texture, I& result);

    std::vector<I> m_isoMipmaps; //< only isotropic mipmaps
    std::vector<std::vector<I>> m_anisoMipmapsWidth; //only anisotropic mipmaps which have reduced width compared to their isotropic counterparts
    std::vector<std::vector<I>> m_anisoMipmapsHeight; //only anisotropic mipmaps which have reduced height compared to their isotropic counterparts

    mipmap_mode_t m_mode;
    bool m_generated;
    bool m_textureSet;
};

template <class I>
Mipmap<I>::Mipmap() :
    m_mode(ISOTROPIC),
    m_generated(false),
    m_textureSet(false)
{}

template <class I>
Mipmap<I>::Mipmap(const I& texture) :
    m_mode(ISOTROPIC),
    m_generated(false),
    m_textureSet(true)
{
    m_isoMipmaps.push_back(texture);
}

template <class I>
void Mipmap<I>::generate(mipmap_mode_t mode, unsigned maxPowReductionLevel)
{
    assert(m_textureSet
           && "Mipmap::generate: no default texture has been set (try using Mipmap::setTexture first)");
    m_generated = true;
    m_mode=mode;
    if(mode == NO_FILTER)
        return;
    //Resizing is done by computing the expected number of mipmaps and comparing it to the max number of mipmaps allowed.

    m_isoMipmaps.resize(   maxPowReductionLevel==0 ? std::floor( std::log2(std::min(m_isoMipmaps[0].width(), m_isoMipmaps[0].height()))+1 )
                                                   : std::min( maxPowReductionLevel+1,
                                                               unsigned(std::floor(std::log2(std::min(  m_isoMipmaps[0].width(),
                                                                                                        m_isoMipmaps[0].height()))))+1 )  );

    for(typename std::vector<I>::iterator it=m_isoMipmaps.begin()+1; it!=m_isoMipmaps.end(); ++it)
    {
        //then successively divide the width and height by 2, computing the average over a 4 pixels window.
        filterDivide2Full(*(it-1), *it);
    }
    //The vector is filled with all isotropic mipmaps.
    //Anisotropic filtering can be enabled by passing MIPMAP_ANISO as the mode parameter of this constructor.
    if(mode==ANISOTROPIC)
    {
        int powReductionLevel; //< sliding max number of reductions by 2, used to cap the reductions
        int w, h, indexIso; //< sliding width, sliding height, index in the isotropic mipmaps vector
        bool cappedNbOfReductions=maxPowReductionLevel==0; //< boolean determining if the number of reductions are capped, for readability

        //First triangular matrix is filled only with mipmaps whose content are mipmaps with successively divided widths compared to isotropic counterparts.
        //For instance, say I have a 256x256 texture (row number 0) with a 64x64 mipmap (row number 2), the columns would be 32x64 (column 0), 16x64 (column 1)...
        //until the max number of reductions is reached or the texture width is 1.

        powReductionLevel=maxPowReductionLevel;
        m_anisoMipmapsWidth.resize(cappedNbOfReductions     ? std::floor( std::log2(m_isoMipmaps[0].height()) )
                                                            : std::min ( maxPowReductionLevel, unsigned(std::floor(std::log2(m_isoMipmaps[0].height()))) ) );
        indexIso=0; //< index used to track the row of the matrix, so we know which isotropic mipmap should be used for the reduction.
        w=m_isoMipmaps[0].width(), h=m_isoMipmaps[0].height();
        //for each row of the matrix, do
        for(typename std::vector<std::vector<I>>::iterator it=m_anisoMipmapsWidth.begin(); it!=m_anisoMipmapsWidth.end(); ++it)
        {
            //pre-compute the size
            (*it).resize(cappedNbOfReductions   ? std::floor( std::log2(w) )
                                                : std::min( --powReductionLevel, int(std::floor(std::log2(w))) ));
            //divide the current mipmap once
            filterDivide2Width(m_isoMipmaps[indexIso++], (*it)[0]);

            //then successively divide the divided mipmaps
            for(typename std::vector<I>::iterator it2=(*it).begin()+1; it2!=(*it).end(); ++it2)
            {
                filterDivide2Width(*(it2-1), *it2);
            }
            //finally, reduce the size of the texture for vector size computations (only w is useful in this loop, but this is for the nice symmetry with the code bellow).
            w/=2;
            h/=2;
        }

        //The second triangular matrix is filled like the first one but the height is divided instead of the width.
        //If the diagonal is set to be the isotropic vector, the concatenation between the first matrix and the transpose of this one is the full mipmap matrix.

        powReductionLevel=maxPowReductionLevel;
        m_anisoMipmapsHeight.resize(cappedNbOfReductions    ? std::floor( std::log2(m_isoMipmaps[0].width()) )
                                                            : std::min( powReductionLevel, int(std::floor(std::log2(m_isoMipmaps[0].width()))) ) );
        indexIso=0;
        w=m_isoMipmaps[0].width(), h=m_isoMipmaps[0].height();
        for(typename std::vector<std::vector<I>>::iterator it=m_anisoMipmapsHeight.begin(); it!=m_anisoMipmapsHeight.end(); ++it)
        {
            (*it).resize(cappedNbOfReductions   ? std::floor( std::log2(h) )
                                                : std::min( --powReductionLevel, int(std::floor(std::log2(h))) ));
            filterDivide2Height(m_isoMipmaps[indexIso++], (*it)[0]);

            for(typename std::vector<I>::iterator it2=(*it).begin()+1; it2!=(*it).end(); ++it2)
            {
                filterDivide2Height(*(it2-1), *it2);
            }
            w/=2;
            h/=2;
        }

    }
}

template <class I>
void Mipmap<I>::revertSetTexture()
{
    if(m_textureSet)
    {
        revertGenerate();
        m_isoMipmaps.erase(m_isoMipmaps.begin());
    }
}

template <class I>
void Mipmap<I>::revertGenerate()
{
    if(m_generated)
    {
        m_isoMipmaps.erase(m_isoMipmaps.begin()+1, m_isoMipmaps.end());
        if(m_mode==ANISOTROPIC)
        {
            m_anisoMipmapsWidth.clear();
            m_anisoMipmapsHeight.clear();
        }
    }
}

template <class I>
void Mipmap<I>::setTexture(const I& texture)
{
    m_generated=false;
    m_textureSet=true;
    if(m_isoMipmaps.size()>0)
        m_isoMipmaps[0]=texture;
    else
        m_isoMipmaps.push_back(texture);
}

template <class I>
size_t Mipmap<I>::numberMipmapsWidth() const
{
    if(m_textureSet)
        return 1;
    if(!m_generated)
        return 0;
    return m_mode == ANISOTROPIC ? m_anisoMipmapsWidth[0].size()+1 : m_isoMipmaps.size();
    //+1 because the first mipmap is the full image, which is stored in m_isoMipmaps[0] but not in m_anisoMipmapsWidth
}

template <class I>
size_t Mipmap<I>::numberMipmapsHeight() const
{
    if(m_textureSet)
        return 1;
    if(!m_generated)
        return 0;
    return m_mode == ANISOTROPIC ? m_anisoMipmapsHeight[0].size()+1 : m_isoMipmaps.size();
}

template <class I>
const I& Mipmap<I>::mipmap(unsigned xPowReduction, unsigned yPowReduction) const
{
    assert(m_mode!=NO_FILTER &&
            "Mipmap::mipmap: cannot call mipmaps with mode set to NO_FILTER (try Mipmap::generate with ISOTROPIC)");
    assert(m_generated &&
            "Mipmap::mipmap: mipmaps have not been generated yet (try using Mipmap::generate)");
    assert((m_mode!=ISOTROPIC || xPowReduction==yPowReduction) &&
            "Mipmap::mipmap: xReduction and yReduction must be identical when using isotropic filtering");
    xPowReduction = std::min(xPowReduction, (unsigned)m_isoMipmaps.size()-1);
    yPowReduction = std::min(yPowReduction, (unsigned)m_isoMipmaps.size()-1);

    if(xPowReduction==yPowReduction)
    {
        return m_isoMipmaps[xPowReduction];
    }

    //The wanted mipmap is in the vector which holds mipmaps with width reduced only if the reduction on x is higher than the one on y ;
    //for access, each line of the anisotropic vector holds images that are already scaled down on x and y.
    //For instance, at line (7, 2), the image is actually scaled down by (7, 9 (7+2)), so we have to take that in account.
    //also, the first column of any anisotropic mipmap vector was made to hold a first reduction (so, that's the -1 in the last index) ;
    //finally, the anisotropic mipmap vector holds lines of different y reductions for the Width one, and vice versa for the Height one, so the indexes may look a little weird.
    return xPowReduction > yPowReduction ? m_anisoMipmapsWidth[yPowReduction][xPowReduction-yPowReduction-1] : m_anisoMipmapsHeight[xPowReduction][yPowReduction-xPowReduction-1];
}

template <class I>
const I& Mipmap<I>::texture() const
{
    assert(m_textureSet &&
            "Mipmap::image(): mipmap must have been given a base texture (use Mipmap::setTexture to give one).");
    return m_isoMipmaps[0];
}

template <class I>
I& Mipmap<I>::texture()
{
    return const_cast<I&>(static_cast<const Mipmap<I>*>(this)->texture());
}

template <class I>
void Mipmap<I>::fullMipmap(I& fullMipmap)
{
    assert(m_textureSet &&
           "Mipmap::fullMipmap: cannot save an empty mipmap (try using Mipmap::setTexture and Mipmap::generate first)");
    //creates a compact image containing all mipmaps computed into fullMipmap.

    //first, determine the width and height of the image, iteratively so we don't bother too much with special cases.
    int w=0, h=0;
    for(typename std::vector<I>::const_iterator it=m_isoMipmaps.begin(); it!=m_isoMipmaps.end(); ++it)
    {
        w+=(*it).width();
        h+=(*it).height();
    }
    fullMipmap.initItk(w, h, true);

    int x0=0, y0=0;
    //then draw the diagonal
    for(typename std::vector<I>::const_iterator it=m_isoMipmaps.begin(); it!=m_isoMipmaps.end(); ++it)
    {
        (*it).for_all_pixels([&] (const typename I::PixelType& pix, int x, int y)
        {
            fullMipmap.pixelAbsolute(x0+x, y0+y) = pix;
        });
        x0+=(*it).width();
        y0+=(*it).height();
    }

    x0=0, y0=0;
    //then draw the rest of the owl in case of anisotropy
    if(m_mode==ANISOTROPIC)
    {
        //superior triangle

        typename std::vector<I>::const_iterator itIso=m_isoMipmaps.begin();
        for(typename std::vector<std::vector<I>>::const_iterator it=m_anisoMipmapsWidth.begin(); it!=m_anisoMipmapsWidth.end(); ++it)
        {
            x0+=(*itIso).width();
            int x1=x0;

            for(typename std::vector<I>::const_iterator it2=(*it).begin(); it2!=(*it).end(); ++it2)
            {
                (*it2).for_all_pixels([&] (const typename I::PixelType& pix, int x, int y)
                {
                    fullMipmap.pixelAbsolute(x1+x, y0+y) = pix;
                });
                x1+=(*it2).width();
            }

            y0+=(*itIso).height();
            ++itIso;
        }

        //inferior triangle

        x0=0, y0=0;
        itIso=m_isoMipmaps.begin();
        for(typename std::vector<std::vector<I>>::const_iterator it=m_anisoMipmapsHeight.begin(); it!=m_anisoMipmapsHeight.end(); ++it)
        {
            y0+=(*itIso).height();
            int y1=y0;

            for(typename std::vector<I>::const_iterator it2=(*it).begin(); it2!=(*it).end(); ++it2)
            {
                (*it2).for_all_pixels([&] (const typename I::PixelType& pix, int x, int y)
                {
                    fullMipmap.pixelAbsolute(x0+x, y1+y) = pix;
                });
                y1+=(*it2).height();
            }

            x0+=(*itIso).width();
            ++itIso;
        }
    }
}

template <class I>
void Mipmap<I>::filterDivide2Width(const I& texture, I& result)
{
    //an average filter, nothing special.
    result.initItk(texture.width()/2, texture.height());
    result.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
    {
        pix=(texture.pixelAbsolute(2*x, y) + texture.pixelAbsolute(2*x + 1, y)) * 0.5;
    });
}

template <class I>
void Mipmap<I>::filterDivide2Height(const I& texture, I& result)
{
    result.initItk(texture.width(), texture.height()/2);
    result.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
    {
        pix=(texture.pixelAbsolute(x, 2*y) + texture.pixelAbsolute(x, 2*y + 1)) * 0.5;
    });
}

template <class I>
void Mipmap<I>::filterDivide2Full(const I& texture, I& result)
{
    result.initItk(texture.width()/2, texture.height()/2);
    result.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
    {
        pix=(texture.pixelAbsolute(2*x, 2*y)
             + texture.pixelAbsolute(2*x + 1, 2*y)
             + texture.pixelAbsolute(2*x, 2*y + 1)
             + texture.pixelAbsolute(2*x + 1, 2*y + 1)) * 0.25;
    });
}

using ImageGrayb = ImageCommon<ImageGrayBase<bool>, false>;
class MipmapBooleanImage : public Mipmap< ImageGrayb >
{
public:

    MipmapBooleanImage();
    MipmapBooleanImage(const ImageGrayb& texture);

protected:

    void filterDivide2Width(const ImageGrayb& texture, ImageGrayb& result);
    void filterDivide2Height(const ImageGrayb& texture, ImageGrayb& result);
    void filterDivide2Full(const ImageGrayb& texture, ImageGrayb& result);
};

MipmapBooleanImage::MipmapBooleanImage():
    Mipmap<ImageGrayb>()
{}

MipmapBooleanImage::MipmapBooleanImage(const ImageGrayb& texture):
    Mipmap<ImageGrayb>(texture)
{}

void MipmapBooleanImage::filterDivide2Width(const ImageGrayb& texture, ImageGrayb& result)
{
    //sets the new cell to true if any of the interpolating cell is true.
    result.initItk(texture.width()/2, texture.height());
    result.for_all_pixels([&] (bool &pix, int x, int y)
    {
        pix = texture.pixelAbsolute(2*x, y) || texture.pixelAbsolute(2*x+1, y);
    });
}

void MipmapBooleanImage::filterDivide2Height(const ImageGrayb& texture, ImageGrayb& result)
{
    result.initItk(texture.width(), texture.height()/2);
    result.for_all_pixels([&] (bool &pix, int x, int y)
    {
        pix = texture.pixelAbsolute(x, 2*y) || texture.pixelAbsolute(x, 2*y+1);
    });
}

void MipmapBooleanImage::filterDivide2Full(const ImageGrayb& texture, ImageGrayb& result)
{
    result.initItk(texture.width()/2, texture.height()/2);
    result.for_all_pixels([&] (bool &pix, int x, int y)
    {
        pix = texture.pixelAbsolute(2*x, 2*y)
            || texture.pixelAbsolute(2*x+1, 2*y)
            || texture.pixelAbsolute(2*x, 2*y+1)
            || texture.pixelAbsolute(2*x+1, 2*y+1);
    });
}

template<typename I>
class MipmapBitmask : public Mipmap<I>
{
public:

    MipmapBitmask();
    MipmapBitmask(const I& texture);

protected:

    void filterDivide2Width(const I& texture, I& result);
    void filterDivide2Height(const I& texture, I& result);
    void filterDivide2Full(const I& texture, I& result);
};

template<typename I>
MipmapBitmask<I>::MipmapBitmask() :
    Mipmap<I>()
{}

template<typename I>
MipmapBitmask<I>::MipmapBitmask(const I& texture) :
    Mipmap<I>(texture)
{}

template<typename I>
void MipmapBitmask<I>::filterDivide2Width(const I& texture, I& result)
{
    result.initItk(texture.width()/2, texture.height());
    result.for_all_pixels([&] (bool &pix, int x, int y)
    {
        pix = texture.pixelAbsolute(2*x, y) | texture.pixelAbsolute(2*x+1, y);
    });
}

template<typename I>
void MipmapBitmask<I>::filterDivide2Height(const I& texture, I& result)
{
    result.initItk(texture.width(), texture.height()/2);
    result.for_all_pixels([&] (bool &pix, int x, int y)
    {
        pix = texture.pixelAbsolute(x, 2*y) | texture.pixelAbsolute(x, 2*y+1);
    });
}

template<typename I>
void MipmapBitmask<I>::filterDivide2Full(const I& texture, I& result)
{
    result.initItk(texture.width()/2, texture.height()/2);
    result.for_all_pixels([&] (bool &pix, int x, int y)
    {
        pix = texture.pixelAbsolute(2*x, 2*y)
            | texture.pixelAbsolute(2*x+1, 2*y)
            | texture.pixelAbsolute(2*x, 2*y+1)
            | texture.pixelAbsolute(2*x+1, 2*y+1);
    });
}

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
    MipmapCEContent(const Mipmap<I> &contentColor, const Mipmap<ImageGrayd> &patchAlpha);

    void generate(mipmap_mode_t mode=ISOTROPIC, unsigned maxPowReductionLevel=0);
};

template<typename I>
MipmapCEContent<I>::MipmapCEContent(const Mipmap<I> &contentColor, const Mipmap<ImageGrayd> &patchAlpha) :
    Mipmap<I>()
{
    assert(contentColor.mode()==patchAlpha.mode() && contentColor.mipmap(0, 0).size()==patchAlpha.mipmap(0, 0).size()
           && contentColor.numberMipmapsWidth() == patchAlpha.numberMipmapsWidth()
           && contentColor.numberMipmapsHeight() == patchAlpha.numberMipmapsHeight()
           && "MipmapCEContent(): contentColor mipmaps and patchAlpha mipmaps don't seem to come from the same process (different mipmap modes or differents sizes)");
    this->m_mode=contentColor.mode();
    this->m_textureSet=true;
    this->m_generated=true;

    unsigned i, j, maxIterations;
    unsigned xMin, yMin, xMax, yMax;
    auto emplaceMipmap = [&] ()
    {
        const I& mipmapContentColor=contentColor.mipmap(i, j);
        const ImageGrayd &mipmapPatchAlpha=patchAlpha.mipmap(i, j);
        I* mipmapContentBox=NULL;

        //slighty awkward assignement, corresponds to mipmap(i, j)
        if(i == j) //isotropic
            mipmapContentBox = &this->m_isoMipmaps[i];
        else
            mipmapContentBox = i > j ? &this->m_anisoMipmapsWidth[j][i-j-1] : &this->m_anisoMipmapsHeight[i][j-i-1];

        //first iteration: find the edges of the box using the patch's alpha.
        //note how patches can bleed from one edge of the image to another due to periodicity.

        xMin=mipmapPatchAlpha.width()-1;
        yMin=mipmapPatchAlpha.height()-1;
        xMax=0;
        yMax=0;
        mipmapPatchAlpha.for_all_pixels([&] (typename ImageGrayd::PixelType &pix, int x, int y)
        {
            if(pix>0) //meaning there's a patch portion
            {
                xMin=std::min(xMin, (unsigned)x);
                yMin=std::min(yMin, (unsigned)y);
                xMax=std::max(xMax, (unsigned)x);
                yMax=std::max(xMax, (unsigned)y);
            }
        });

        //we have found a non periodic bounding box.
        //Now we need to check if the bounding box is aperiodic in width.
        //TODO: if a bug is found, it could come from this part. It's easy to make a mistake.
        //TODO: this part of the code WILL also be problematic if some periodic boxes take the entire width.
        if(xMin == 0 && xMax == mipmapPatchAlpha.width()-1)
        { //this block occurs when the bounding box is suspected to be periodic in width.
            //check discontinuities in image width (an empty column)
            int oldXMax = xMax;
            for(int x=xMin; x<=xMax; ++x)
            {
                int y;
                for(y=yMin; y<=yMax && mipmapPatchAlpha.pixelAbsolute(x, y)==0; ++y);
                if(y>yMax) //discontinuity found
                    xMax=x-1; //assign xMax to the last full row we found, breaking the loop
            }
            if(oldXMax != xMax) //if false, it means the content was simply taking the entire width
            {
                //iterate the other way until a discontinuity is found
                for(int x=oldXMax; x>=xMin; --x)
                {
                    int y;
                    for(y=yMin; y<=yMax && mipmapPatchAlpha.pixelAbsolute(x, y)==0; ++y);
                    if(y>yMax) //discontinuity found
                        xMin=x+1; //assign xMin to the last full row we found, breaking the loop
                }
            }
        }
        //from here onwards, the bounding box shall take periodicity in width in account.
        //now for the periodicity in height, it's the same iteration, but in height, with y instead of x, and a modulo on x.
        if(yMin == 0 && yMax == mipmapPatchAlpha.height()-1)
        {
            int oldYMax = yMax;
            //notice the new iterator: we use != to xMax+1 and the iteration is periodic on the width.
            for(int y=yMin; y<=yMax; ++y)
            {
                int x;
                for(x=xMin; x!=xMax+1 && mipmapPatchAlpha.pixelAbsolute(x, y)==0; x = (x+1)%mipmapPatchAlpha.width());
                if(x==xMax+1) //discontinuity found (don't use >)
                    yMax=y-1;
            }
            if(oldYMax != yMax)
            {
                for(int y=oldYMax; y>=yMin; --y)
                {
                    int x;
                    for(x=xMin; x!=xMax+1 && mipmapPatchAlpha.pixelAbsolute(x, y)==0; x = (x+1)%mipmapPatchAlpha.width());
                    if(x==xMax+1)
                        yMin=y+1;
                }
            }
        }
        //from here onwards we have found the correctly periodic bounding box.
        //TODO: this code should be copied, or even exported, to the place we compute the boxes for the patches.
        //Then, we could simply pass every boxes we computed to this function.

        //Allocate the new mipmap

        mipmapContentBox->initItk(xMin<xMax ? xMax-xMin+1 : mipmapPatchAlpha.width()-xMin + xMax,
                                  yMin<yMax ? yMax-yMin+1 : mipmapPatchAlpha.height()-yMin + yMax, true);
        //Fill up the new mipmap
        mipmapContentBox->for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
        {
            pix = mipmapContentColor.pixelAbsolute((xMin+x)%mipmapContentColor.width(),
                                                   (yMin+y)%mipmapContentColor.height());
        });
    };

    //allocation

    maxIterations = std::max(contentColor.numberMipmapsWidth(), contentColor.numberMipmapsHeight());
    this->m_isoMipmaps.resize(maxIterations);
    if(contentColor.mode() == Mipmap<I>::ANISOTROPIC)
    {
        this->m_anisoMipmapsWidth.resize(contentColor.numberMipmapsWidth()-1);
        this->m_anisoMipmapsHeight.resize(contentColor.numberMipmapsHeight()-1);
    }

    if(contentColor.mode()==Mipmap<I>::ISOTROPIC)
    for(i=0, j=0; i<maxIterations; ++i, ++j)
        emplaceMipmap();

    else if(contentColor.mode()==Mipmap<I>::ANISOTROPIC)
        for(i=0; i<contentColor.numberMipmapsWidth(); ++i)
            for(j=0; j<contentColor.numberMipmapsHeight(); ++i)
                emplaceMipmap();


}

template<typename I>
void MipmapCEContent<I>::generate(mipmap_mode_t mode, unsigned maxPowReductionLevel)
{
    (void) mode;
    (void) maxPowReductionLevel;
    return;
}

#endif //__MIPMAP__
