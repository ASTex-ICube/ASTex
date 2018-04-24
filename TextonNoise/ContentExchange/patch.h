#ifndef __PATCH__H__
#define __PATCH__H__

#include <ASTex/image_gray.h>
#include "mipmap.h"
//#include "content.h"

namespace ASTex
{

namespace ContentExchange
{

template<typename I>
class Content;

using PixelPos = itk::Index<2>;
using ImageAlphad = ImageGrayd;

//Mipmap

class MipmapCEPatch : public Mipmap<ImageAlphad>
{
public:

    typedef itk::Index<2> PixelPos;

    MipmapCEPatch();
    MipmapCEPatch(const ImageAlphad &patchAlpha);

    void generate();

    PixelPos originAt(unsigned xPowReduction, unsigned yPowReduction) const;

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

void MipmapCEPatch::generate()
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
    IO::save01_in_u8(fullMipmapAlphaImage, "fullMipmapAlphaImage.png");

    //origins vector: it's okay for it to be redundant in case of isotropy since it's not heavy
    m_pixelOriginMap.resize(this->numberMipmapsWidth());
    for(std::vector<std::vector<PixelPos>>::iterator it=m_pixelOriginMap.begin(); it!=m_pixelOriginMap.end(); ++it)
    {
        (*it).resize(this->numberMipmapsHeight());
    }

    unsigned i, j;
    int xMin, yMin, xMax, yMax;
    int oldXMin, oldXMax, oldYMax;

    static int fffff=0;

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
                std::cout << xMax << std::endl;
                for(x=oldXMin; x<=oldXMax && mipmapAlpha.pixelAbsolute(x, y)==0; ++x);
                if(x>oldXMax) //discontinuity found (don't use >)
                    yMax=y-1;
            }
            if(oldYMax != yMax)
            {
                for(int y=oldYMax; y>=yMin; --y)
                {
                    int x;
                    for(x=oldXMin; x<oldXMax && mipmapAlpha.pixelAbsolute(x, y)==0; ++x);
                    if(x>oldXMax)
                        yMin=y+1;
                }
            }
        }


        //from here onwards we have found the correctly periodic bounding box.
        //Allocate the new mipmap

        mipmapPatchAlpha.initItk(   xMin<xMax ? xMax-xMin+1 : mipmapAlpha.width()-xMin + xMax,
                                    yMin<yMax ? yMax-yMin+1 : mipmapAlpha.height()-yMin + yMax, true);
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


#include "ASTex/easy_io.h"
        IO::save01_in_u8(mipmapAlpha, "mipmapPatch_" + std::to_string(fffff) + "_number_" + std::to_string(i) + "_" + std::to_string(j) + ".png");
    };

    if(m_mode==NO_FILTER || m_mode==ISOTROPIC)
        for(i=0, j=0; i<m_isoMipmaps.size(); ++i, ++j)
            emplaceMipmap();
    else if(m_mode==ANISOTROPIC)
        for(i=0; i<numberMipmapsWidth(); ++i)
            for(j=0; j<numberMipmapsHeight(); ++j)
                emplaceMipmap();

    ++fffff;

    return;
}

MipmapCEPatch::PixelPos MipmapCEPatch::originAt(unsigned xPowReduction, unsigned yPowReduction) const
{
    size_t xBoundedIndex, yBoundedIndex;
    xBoundedIndex=std::min((size_t)xPowReduction, numberMipmapsWidth()-1);
    yBoundedIndex=std::min((size_t)yPowReduction, numberMipmapsHeight()-1);
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

    const ImageAlphad& mipmap(unsigned i, unsigned j) const {return m_alphaMipmap.mipmap(i, j);}
    const MipmapCEPatch& alphaMipmap() const {return m_alphaMipmap;}
    const Content<I>& contentAt(size_t index) const {return m_alternativeContents[index];}
    Content<I>& contentAt(size_t index) {return m_alternativeContents[index];}
    size_t numberContents() const {return m_alternativeContents.size();}

    //set

    void setAlphaMap(const ImageAlphad& alphaMap, mipmap_mode_t mipmapMode=NO_FILTER);
    void addContent(const Content<I> &content);

private:

    PixelPos m_origin;
    MipmapCEPatch m_alphaMipmap;

    std::vector<Content<I>> m_alternativeContents;
    //and a probability map for choosing them.
};

template<typename I>
Patch<I>::Patch():
    m_origin(),
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
