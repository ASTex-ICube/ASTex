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
    Mipmap(const Mipmap& other);

    //get

    mipmap_mode_t mode() const {return m_mode;}
    unsigned maxPowReductionLevel() const {return m_maxPowReductionLevel;}

    size_t numberMipmapsWidth() const;
    size_t numberMipmapsHeight() const;
    const I& mipmap(unsigned xPowReduction, unsigned yPowReduction) const;
    I& mipmap(unsigned xPowReduction, unsigned yPowReduction);

    const I& texture() const;
    I& texture();

    bool isGenerated() const {return m_generated;}
    bool isTextureSet() const {return m_textureSet;}

    //set

    void setTexture(const I& texture);
    void setMode(mipmap_mode_t mode);
    void setMaxPowReductionLevel(unsigned maxPowReductionLevel=0);

    //misc

    virtual void generate();
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
    unsigned m_maxPowReductionLevel;
    bool m_generated;
    bool m_textureSet;
};

template <class I>
Mipmap<I>::Mipmap() :
    m_mode(ISOTROPIC),
    m_maxPowReductionLevel(0),
    m_generated(false),
    m_textureSet(false)
{}

template <class I>
Mipmap<I>::Mipmap(const I& texture) :
    m_mode(ISOTROPIC),
    m_maxPowReductionLevel(0),
    m_generated(false),
    m_textureSet(true)
{
    this->setTexture(texture);
}

template <class I>
Mipmap<I>::Mipmap(const Mipmap& other):
    m_isoMipmaps(),
    m_anisoMipmapsWidth(),
    m_anisoMipmapsHeight(),
    m_mode(other.mode()),
    m_maxPowReductionLevel(other.maxPowReductionLevel()),
    m_generated(other.isGenerated()),
    m_textureSet(other.isTextureSet())
{
    if(m_textureSet)
    {
        m_isoMipmaps.resize(std::max(other.numberMipmapsWidth(), other.numberMipmapsHeight()));
        m_anisoMipmapsWidth.resize(other.numberMipmapsWidth()-1); //counter-intuitive, but trust me
        m_anisoMipmapsHeight.resize(other.numberMipmapsHeight()-1);
        for(size_t i=0; i<m_anisoMipmapsWidth.size(); ++i)
        {
            m_anisoMipmapsWidth[i].resize(other.numberMipmapsWidth()-1-i);
        }
        for(size_t j=0; j<m_anisoMipmapsHeight.size(); ++j)
        {
            m_anisoMipmapsHeight[j].resize(other.numberMipmapsHeight()-1-j);
        }
        for(size_t i=0; i<other.numberMipmapsWidth(); ++i)
            for(size_t j=0; j<other.numberMipmapsHeight(); ++j)
            {
                I& m=this->mipmap(i, j);
                const I& mOther=other.mipmap(i, j);
                m.initItk(mOther.width(), mOther.height(), false);
                m.copy_pixels(mOther);
            }
    }
}

template <class I>
void Mipmap<I>::generate()
{
    assert(m_textureSet
           && "Mipmap::generate: no default texture has been set (try using Mipmap::setTexture first)");
    m_generated = true;
    if(m_mode == NO_FILTER)
        return;
    //Resizing is done by computing the expected number of mipmaps and comparing it to the max number of mipmaps allowed.

    m_isoMipmaps.resize(   m_maxPowReductionLevel==0 ? std::floor( std::log2(std::min(m_isoMipmaps[0].width(), m_isoMipmaps[0].height()))+1 )
                                                   : std::min( m_maxPowReductionLevel+1,
                                                               unsigned(std::floor(std::log2(std::min(  m_isoMipmaps[0].width(),
                                                                                                        m_isoMipmaps[0].height()))))+1 )  );

    for(typename std::vector<I>::iterator it=m_isoMipmaps.begin()+1; it!=m_isoMipmaps.end(); ++it)
    {
        //then successively divide the width and height by 2, computing the average over a 4 pixels window.
        filterDivide2Full(*(it-1), *it);
    }
    //The vector is filled with all isotropic mipmaps.
    //Anisotropic filtering can be enabled by passing MIPMAP_ANISO as the mode parameter of this constructor.
    if(m_mode==ANISOTROPIC)
    {
        int powReductionLevel; //< sliding max number of reductions by 2, used to cap the reductions
        int w, h, indexIso; //< sliding width, sliding height, index in the isotropic mipmaps vector
        bool cappedNbOfReductions=m_maxPowReductionLevel==0; //< boolean determining if the number of reductions are capped, for readability

        //First triangular matrix is filled only with mipmaps whose content are mipmaps with successively divided widths compared to isotropic counterparts.
        //For instance, say I have a 256x256 texture (row number 0) with a 64x64 mipmap (row number 2), the columns would be 32x64 (column 0), 16x64 (column 1)...
        //until the max number of reductions is reached or the texture width is 1.

        powReductionLevel=m_maxPowReductionLevel;
        m_anisoMipmapsWidth.resize(cappedNbOfReductions     ? std::floor( std::log2(m_isoMipmaps[0].height()) )
                                                            : std::min ( m_maxPowReductionLevel, unsigned(std::floor(std::log2(m_isoMipmaps[0].height()))) ) );
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

        powReductionLevel=m_maxPowReductionLevel;
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
        m_generated=false;
    }
}

template <class I>
void Mipmap<I>::setMode(mipmap_mode_t mode)
{
    if(m_mode!=mode)
        revertGenerate();
    m_mode=mode;
}

template <class I>
void Mipmap<I>::setMaxPowReductionLevel(unsigned maxPowReductionLevel)
{
    if(maxPowReductionLevel!=m_maxPowReductionLevel)
        revertGenerate();
    m_maxPowReductionLevel=maxPowReductionLevel;
}

template <class I>
void Mipmap<I>::setTexture(const I& texture)
{
    m_generated=false;
    //revertGenerate(); //Yes? No?
    m_textureSet=true;
    if(m_isoMipmaps.size()>0)
    {
        //because apparently = and the copy constructor can't be bothered to deep copy half the time
        m_isoMipmaps[0].initItk(texture.width(), texture.height());
        m_isoMipmaps[0].copy_pixels(texture);
    }
    else
    {
        //same here
        I idiotTexture;
        idiotTexture.initItk(texture.width(), texture.height());
        idiotTexture.copy_pixels(texture);
        m_isoMipmaps.push_back(idiotTexture);
    }
}

template <class I>
size_t Mipmap<I>::numberMipmapsWidth() const
{
    if(!m_generated && m_textureSet)
        return 1;
    if(!m_generated)
        return 0;
    return m_mode == ANISOTROPIC ? m_anisoMipmapsWidth[0].size()+1 : m_isoMipmaps.size();
    //+1 because the first mipmap is the full image, which is stored in m_isoMipmaps[0] but not in m_anisoMipmapsWidth
}

template <class I>
size_t Mipmap<I>::numberMipmapsHeight() const
{
    if(!m_generated && m_textureSet)
        return 1;
    if(!m_generated)
        return 0;
    return m_mode == ANISOTROPIC ? m_anisoMipmapsHeight[0].size()+1 : m_isoMipmaps.size();
}

template <class I>
const I& Mipmap<I>::mipmap(unsigned xPowReduction, unsigned yPowReduction) const
{
    assert((m_mode!=NO_FILTER || (xPowReduction==0 && yPowReduction==0)) &&
            "Mipmap::mipmap: cannot use mipmaps with mode set to NO_FILTER (try Mipmap::setMode with ISOTROPIC/ANISOTROPIC)");
    assert(m_generated &&
            "Mipmap::mipmap: mipmaps have not been generated yet (try using Mipmap::generate)");
    assert((m_mode!=ISOTROPIC || xPowReduction==yPowReduction) &&
            "Mipmap::mipmap: xReduction and yReduction must be identical when using isotropic filtering (or try Mipmap::setMode with ANISOTROPIC)");
    xPowReduction = std::min(xPowReduction, (unsigned)m_isoMipmaps.size()-1);
    yPowReduction = std::min(yPowReduction, (unsigned)m_isoMipmaps.size()-1);

    if(xPowReduction==yPowReduction)
        return m_isoMipmaps[xPowReduction];

    //The wanted mipmap is in the vector which holds mipmaps with width reduced only if the reduction on x is higher than the one on y ;
    //for access, each line of the anisotropic vector holds images that are already scaled down on x and y.
    //For instance, at line (7, 2), the image is actually scaled down by (7, 9 (7+2)), so we have to take that in account.
    //also, the first column of any anisotropic mipmap vector was made to hold a first reduction (so, that's the -1 in the last index) ;
    //finally, the anisotropic mipmap vector holds lines of different y reductions for the Width one, and vice versa for the Height one, so the indexes may look a little weird.
    return xPowReduction > yPowReduction
            ? m_anisoMipmapsWidth[yPowReduction][xPowReduction-yPowReduction-1]
            : m_anisoMipmapsHeight[xPowReduction][yPowReduction-xPowReduction-1];
}

template <class I>
I& Mipmap<I>::mipmap(unsigned xPowReduction, unsigned yPowReduction)
{
    return const_cast<I&>(static_cast<const Mipmap<I>*>(this)->mipmap(xPowReduction, yPowReduction));
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
    result.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
    {
        pix = texture.pixelAbsolute(2*x, y) | texture.pixelAbsolute(2*x+1, y);
    });
}

template<typename I>
void MipmapBitmask<I>::filterDivide2Height(const I& texture, I& result)
{
    result.initItk(texture.width(), texture.height()/2);
    result.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
    {
        pix = texture.pixelAbsolute(x, 2*y) | texture.pixelAbsolute(x, 2*y+1);
    });
}

template<typename I>
void MipmapBitmask<I>::filterDivide2Full(const I& texture, I& result)
{
    result.initItk(texture.width()/2, texture.height()/2);
    result.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
    {
        pix = texture.pixelAbsolute(2*x, 2*y)
            | texture.pixelAbsolute(2*x+1, 2*y)
            | texture.pixelAbsolute(2*x, 2*y+1)
            | texture.pixelAbsolute(2*x+1, 2*y+1);
    });
}

#endif //__MIPMAP__
