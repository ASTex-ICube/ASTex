#ifndef __MIPMAP__H__
#define __MIPMAP__H__

#include <cmath>
#include <ASTex/special_io.h>
#include <ASTex/utils.h>
#include <ASTex/colorspace_filters.h>
#include <ASTex/mask.h>
#include <type_traits>

using namespace ASTex;

template <class I>
class Mipmap
{
public:

    typedef enum {MIPMAP_ISO=0, MIPMAP_ANISO=1} mipmap_mode_t;

    Mipmap(const I& texture, mipmap_mode_t mode=MIPMAP_ISO, int maxReductionLevel=0);

    const I& mipmap(int xPowReduction, int yPowReduction);

    void fullMipmap(I& fullMipmap);

private:

    void filterDivide2Width(const I& texture, I& result);
    void filterDivide2Height(const I& texture, I& result);
    void filterDivide2Full(const I& texture, I& result);

    std::vector<I> m_isoMipmaps; //< only isotropic mipmaps
    std::vector<std::vector<I>> m_anisoMipmapsWidth; //only anisotropic mipmaps which have reduced width compared to their isotropic counterparts
    std::vector<std::vector<I>> m_anisoMipmapsHeight; //only anisotropic mipmaps which have reduced height compared to their isotropic counterparts

    mipmap_mode_t m_mode;
};

template <class I>
Mipmap<I>::Mipmap(const I& texture, mipmap_mode_t mode, int maxPowReductionLevel) :
    m_mode(mode)
{
    //Resizing is done by computing the expected number of mipmaps and comparing it to the max number of mipmaps allowed.

    m_isoMipmaps.resize(   maxPowReductionLevel<=0 ? std::floor( std::log2(std::min(texture.width(), texture.height()))+1 )
                                                   : std::min( maxPowReductionLevel+1, int(std::floor(std::log2(std::min(texture.width(), texture.height()))))+1 )  );

    //Set the first mipmap to the given image
    m_isoMipmaps[0].initItk(texture.width(), texture.height());
    m_isoMipmaps[0].copy_pixels(texture);
    for(typename std::vector<I>::iterator it=m_isoMipmaps.begin()+1; it!=m_isoMipmaps.end(); ++it)
    {
        //then successively divide the width and height by 2, computing the average over a 4 pixels window.
        filterDivide2Full(*(it-1), *it);
    }
    //The vector is filled with all isotropic mipmaps.
    //Anisotropic filtering can be enabled by passing MIPMAP_ANISO as the mode parameter of this constructor.
    if(mode==MIPMAP_ANISO)
    {
        int powReductionLevel; //< sliding max number of reductions by 2, used to cap the reductions
        int w, h, indexIso; //< sliding width, sliding height, index in the isotropic mipmaps vector
        bool cappedNbOfReductions=maxPowReductionLevel<=0; //< boolean determining if the number of reductions are capped, for readability

        //First triangular matrix is filled only with mipmaps whose content are mipmaps with successively divided widths compared to isotropic counterparts.
        //For instance, say I have a 256x256 texture (row number 0) with a 64x64 mipmap (row number 2), the columns would be 32x64 (column 0), 16x64 (column 1)...
        //until the max number of reductions is reached or the texture width is 1.

        powReductionLevel=maxPowReductionLevel;
        m_anisoMipmapsWidth.resize(cappedNbOfReductions     ? std::floor( std::log2(texture.height()) )
                                                            : std::min ( maxPowReductionLevel, int(std::floor(std::log2(texture.height()))) ) );
        indexIso=0; //< index used to track the row of the matrix, so we know which isotropic mipmap should be used for the reduction.
        w=texture.width(), h=texture.height();
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
        m_anisoMipmapsHeight.resize(cappedNbOfReductions    ? std::floor( std::log2(texture.width()) )
                                                            : std::min( powReductionLevel, int(std::floor(std::log2(texture.width()))) ) );
        indexIso=0;
        w=texture.width(), h=texture.height();
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
const I& Mipmap<I>::mipmap(int xPowReduction, int yPowReduction)
{
    assert(xPowReduction>=0 && yPowReduction>=0 && "Mipmap::mipmap: xReduction and yReduction must be >= 0");
    assert((m_mode!=MIPMAP_ISO || xPowReduction==yPowReduction) && "Mipmap::mipmap: xReduction and yReduction must be identical when using isotropic filtering");
    xPowReduction = std::min(xPowReduction, (int)m_isoMipmaps.size()-1);
    yPowReduction = std::min(yPowReduction, (int)m_isoMipmaps.size()-1);

    if(xPowReduction==yPowReduction)
    {
        return m_isoMipmaps[xPowReduction];
    }

    //The wanted mipmap is in the vector which holds mipmaps with width reduced only if the reduction on x is higher than the one on y ;
    //for access, each line of the anisotropic vector holds images that are already scaled down on x and y.
    //For instance, at line (7, 2), the image is actually scaled down by (7, 9), so we have to take that in account.
    //also, the first column of any anisotropic mipmap vector was made to hold a first reduction (so, that's the -1 in the last index) ;
    //finally, the anisotropic mipmap vector holds lines of different y reductions for the Width one, and vice versa for the Height one, so the indexes may look a little weird.
    return xPowReduction > yPowReduction ? m_anisoMipmapsWidth[yPowReduction][xPowReduction-yPowReduction-1] : m_anisoMipmapsHeight[xPowReduction][yPowReduction-xPowReduction-1];
}

template <class I>
void Mipmap<I>::filterDivide2Width(const I& texture, I& result)
{
    //an average filter, nothing special.
    result.initItk(texture.width()/2, texture.height());
    result.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
    {
        pix=(texture.pixelAbsolute(2*x, y) + texture.pixelAbsolute(2*x + 1, y)) * (1.0/2.0);
    });
}

template <class I>
void Mipmap<I>::filterDivide2Height(const I& texture, I& result)
{
    result.initItk(texture.width(), texture.height()/2);
    result.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
    {
        pix=(texture.pixelAbsolute(x, 2*y) + texture.pixelAbsolute(x, 2*y + 1)) * (1.0/2.0);
    });
}

template <class I>
void Mipmap<I>::filterDivide2Full(const I& texture, I& result)
{
    result.initItk(texture.width()/2, texture.height()/2);
    result.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
    {
        pix=(texture.pixelAbsolute(2*x, 2*y) + texture.pixelAbsolute(2*x + 1, 2*y) + texture.pixelAbsolute(2*x, 2*y + 1) + texture.pixelAbsolute(2*x + 1, 2*y + 1)) * (1.0/4.0); //TODO: /4.0 doesn't work if I is ImageRGBd, but *(1.0/4.0) does.
    });
}

template <class I>
void Mipmap<I>::fullMipmap(I& fullMipmap)
{
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
    if(m_mode==MIPMAP_ANISO)
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

#endif //__MIPMAP__
