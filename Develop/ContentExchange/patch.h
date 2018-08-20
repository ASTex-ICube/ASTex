#ifndef __CTEXCH_PATCH__H__
#define __CTEXCH_PATCH__H__

#include <ASTex/image_gray.h>
#include "mipmap.h"
//#include "content.h"

namespace ASTex
{

namespace ContentExchange
{

//Contains MipmapCEPatch and Patch classes.

template<typename I>
class Content;

using PixelPos = itk::Index<2>;
using ImageAlphad = ImageGrayd;

//Mipmap

/**
 * @brief The MipmapCEPatch class is a derivation of the mipmap class for alpha,
 * where origins of each patches in the input textures are computed.
 * The class expects an "alpha" map, where the pixel is set to 1 if it corresponds to a certain patch, 0 otherwise.
 * One should mostly use the above Patch class to access it.
 */
class MipmapCEPatch : public Mipmap<ImageAlphad>
{
public:

    typedef itk::Index<2> PixelPos;

    /**
     * @brief MipmapCEPatch if using this constructor, use setTexture(ImageAlphad) to set a texture.
     */
    MipmapCEPatch();
    MipmapCEPatch(const ImageAlphad &patchAlpha);

    /**
     * @brief generate bakes the mipmap with settings given by setTexture(), setMode() and setMaxPowReductionLevel().
     * @param periodicity determines if patch bounding boxes is computed across periodic boundaries or not.
     * This is a critical setting, because it could affect blending as well as what can be pre-computed for the GPU,
     * since displacement computations are easier in a non-periodic environment, but blending computations aren't.
     */
    void generate(bool periodicity = true);

    /**
     * @brief originAt
     * @pre generate() MUST have been called before calling this function.
     * @param xPowReduction specifies the level of reduction on width.
     * @param yPowReduction specifies the level of reduction on height.
     * @return the origin of (this) patch, for the given level of reduction.
     */
    PixelPos originAt(unsigned k, unsigned l) const;

private:

    std::vector<std::vector<PixelPos>> m_pixelOriginMap;
};

MipmapCEPatch::MipmapCEPatch():
    Mipmap<ImageAlphad>(), //sets m_textureSet to true
    m_pixelOriginMap()
{}

MipmapCEPatch::MipmapCEPatch(const ImageAlphad &patchAlpha):
    Mipmap<ImageAlphad>(patchAlpha), //sets m_textureSet to true
    m_pixelOriginMap()
{}

void MipmapCEPatch::generate(bool periodicity)
{
    if(!isTextureSet())
    {
        std::cerr << "Warning: MipmapCEPatch::generate: texture must be set before generating" << std::endl;
        return;
    }
    Mipmap<ImageAlphad>::generate(); //will set m_generated to true and generate a big alpha mipmap map
    //debug
    ImageAlphad fullMipmapAlphaImage;
    this->fullMipmap(fullMipmapAlphaImage);

    //origins vector: it's okay for it to be redundant in case of isotropy since it's not heavy
    m_pixelOriginMap.resize(this->numberMipmapsWidth());
    for(std::vector<std::vector<PixelPos>>::iterator it=m_pixelOriginMap.begin(); it!=m_pixelOriginMap.end(); ++it)
    {
        (*it).resize(this->numberMipmapsHeight());
    }

    unsigned i, j;
    int xMin, yMin, xMax, yMax;
    int oldXMin, oldXMax, oldYMax;

    auto emplaceMipmap = [&] ()
    {
        ImageAlphad &mipmapAlpha=mipmap(i, j); //the mipmap we are going to change in the end
        ImageAlphad mipmapPatchAlpha; //the mipmap we are going to compute

        //first iteration: find the edges of the box using the patch's alpha.
        //note how patches can bleed from one edge of the image to another due to periodicity.

        xMin=mipmapAlpha.width()-1;
        yMin=mipmapAlpha.height()-1;
        xMax=0;
        yMax=0;
        mipmapAlpha.for_all_pixels([&] (const typename ImageAlphad::PixelType &pix, int x, int y)
        {
            if(pix>0) //meaning there's a patch portion
            {
                xMin=std::min(xMin, x);
                yMin=std::min(yMin, y);
                xMax=std::max(xMax, x);
                yMax=std::max(yMax, y);
            }
        });

        oldXMin=xMin;
        oldXMax=xMax;
        oldYMax=yMax;

        //we have found a non periodic bounding box.
        //Now we need to check if the bounding box is aperiodic in width.
        //TODO: if a bug is found, it could come from this part. It's easy to make a mistake.
        //TODO: this part of the code WILL also be problematic if some periodic boxes take the entire width.
        if(periodicity)
        {
            if(xMin == 0 && xMax == mipmapAlpha.width()-1)
            { //this block occurs when the bounding box is suspected to be periodic in width.
                //check discontinuities in image width (an empty column)
                oldXMax = xMax;
                for(int x=xMin; x<=xMax; ++x)
                {
                    int y;
                    for(y=yMin; y<=yMax && mipmapAlpha.pixelAbsolute(x, y)==0; ++y);
                    if(y>yMax) //discontinuity found
                        xMax=x-1; //assign xMax to the last full row we found, breaking the loop
                }
                if(oldXMax != xMax) //if false, it means the content was simply taking the entire width
                {
                    //iterate the other way until a discontinuity is found
                    for(int x=oldXMax; x>=xMin; --x)
                    {
                        int y;
                        for(y=yMin; y<=yMax && mipmapAlpha.pixelAbsolute(x, y)==0; ++y);
                        if(y>yMax) //discontinuity found
                            xMin=x+1; //assign xMin to the last full row we found, breaking the loop
                    }
                }
            }
            //from here onwards, the bounding box shall take periodicity in width in account.
            //now for the periodicity in height, it's the same iteration, but in height, with y instead of x, and a modulo on x.
            if(yMin == 0 && yMax == mipmapAlpha.height()-1)
            {
                int oldYMax = yMax;
                //notice the new iterator: we use != to xMax+1 and the iteration is periodic on the width.
                for(int y=yMin; y<=yMax; ++y)
                {
                    int x;
                    for(x=oldXMin; x<=oldXMax && mipmapAlpha.pixelAbsolute(x, y)==0; ++x);
                    if(x>oldXMax) //discontinuity found (don't use >)
                        yMax=y-1;
                }
                if(oldYMax != yMax)
                {
                    for(int y=oldYMax; y>=yMin; --y)
                    {
                        int x;
                        for(x=oldXMin; x<=oldXMax && mipmapAlpha.pixelAbsolute(x, y)==0; ++x);
                        if(x>oldXMax)
                            yMin=y+1;
                    }
                }
            }
        }


        //from here onwards we have found the correctly periodic bounding box.
        //Allocate the new mipmap

        mipmapPatchAlpha.initItk(   xMin<xMax ? xMax-xMin+1 : mipmapAlpha.width()-xMin + xMax+1,
                                    yMin<yMax ? yMax-yMin+1 : mipmapAlpha.height()-yMin + yMax+1, true);
        //Fill up the new mipmap

        mipmapPatchAlpha.for_all_pixels([&] (ImageAlphad::PixelType &pix, int x, int y)
        {
            pix=mipmapAlpha.pixelAbsolute((xMin+x)%mipmapAlpha.width(), (yMin+y)%mipmapAlpha.height());
        });
        PixelPos origin;
        origin[0]=xMin;
        origin[1]=yMin;
        m_pixelOriginMap[i][j]=origin;
        //emplace the new mipmap
        mipmapAlpha.initItk(mipmapPatchAlpha.width(), mipmapPatchAlpha.height());
        mipmapAlpha.copy_pixels(mipmapPatchAlpha);
    };

    if(m_mode==NO_FILTER || m_mode==ISOTROPIC)
        for(i=0, j=0; i<m_isoMipmaps.size(); ++i, ++j)
            emplaceMipmap();
    else if(m_mode==ANISOTROPIC)
        for(i=0; i<numberMipmapsWidth(); ++i)
            for(j=0; j<numberMipmapsHeight(); ++j)
                emplaceMipmap();

    return;
}

MipmapCEPatch::PixelPos MipmapCEPatch::originAt(unsigned k, unsigned l) const
{
    size_t xBoundedIndex, yBoundedIndex;
    xBoundedIndex=std::min((size_t)k, numberMipmapsWidth()-1);
    yBoundedIndex=std::min((size_t)l, numberMipmapsHeight()-1);
    return m_pixelOriginMap[xBoundedIndex][yBoundedIndex];
}

//Patch

template<typename I>
class Patch
{
public:

    Patch();

    //get

    using ImageAlphad = ImageGrayd;

    /**
     * @brief mipmap
     * @param k
     * @param l
     * @return the alpha mipmap of reduction level k, l
     */
    const ImageAlphad& mipmap(unsigned k, unsigned l) const {return m_alphaMipmap.mipmap(k, l);}

    /**
     * @brief alphaMap
     * @param k
     * @param l
     * @return the alpha mipmap of reduction level 0, 0
     */
    const ImageAlphad& alphaMap() const {return m_alphaMipmap.texture();}

    /**
     * @brief alphaMipmap
     * @return the alpha mipmap for this patch.
     */
    const MipmapCEPatch& alphaMipmap() const {return m_alphaMipmap;}

    /**
     * @brief contentAt
     * @param index
     * @return returns content number index.
     * @pre contents must have been set by the parent PatchProcessor, or it will produce undefined behavior.
     */
    const Content<I>& contentAt(size_t index) const {return m_alternativeContents[index];}
    Content<I>& contentAt(size_t index) {return m_alternativeContents[index];}

    /**
     * @brief originAt
     * @param k
     * @param l
     * @return alphaMipmap().originAt(k, l), which corresponds to the origin of the patch for reduction level (k, l).
     */
    PixelPos originAt(unsigned k, unsigned l) const {return m_alphaMipmap.originAt(k, l);}

    /**
     * @brief numberContents
     * @return the number of contents. This number bounds contentAt().
     */
    size_t nbContents() const {return m_alternativeContents.size();}

    //set

    /**
     * @brief setAlphaMap patches are defined by their position in the input.
     * The alphaMap structure corresponds to a texture entity of the same size as the input, where :
     * alphaMap(x) = 0 if x doesn't belong to this patch, and 1 if it does.
     *
     * This function builds the entire mipmap required with mode mipmapMode.
     * @param alphaMap
     * @param mipmapMode
     */
    void setAlphaMap(const ImageAlphad& alphaMap, mipmap_mode_t mipmapMode=NO_FILTER);

    /**
     * @brief addContent appends a content to the content collection.
     * Since ContentExchange was build to store contents, this could be any image,
     * as long as it fits inside the patch.
     * @param content
     */
    void addContent(const Content<I> &content);

private:

    MipmapCEPatch m_alphaMipmap;

    std::vector<Content<I>> m_alternativeContents;
    //and a probability map for choosing them.
};

template<typename I>
Patch<I>::Patch():
    m_alphaMipmap(),
    m_alternativeContents()
{}

template<typename I>
void Patch<I>::setAlphaMap(const ImageAlphad& alphaMap, mipmap_mode_t mipmapMode)
{
    m_alphaMipmap.setTexture(alphaMap);
    m_alphaMipmap.setMode(mipmapMode);
    m_alphaMipmap.generate();
}

template<typename I>
void Patch<I>::addContent(const Content<I> &content)
{
    m_alternativeContents.push_back(content);
}

} //contentexchange namespace

}


#endif
