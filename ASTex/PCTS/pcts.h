#ifndef __PCTS_H__
#define __PCTS_H__

#include "ASTex/mipmap.h"
#include "ASTex/easy_io.h"

namespace ASTex
{

using ImageIndex2 = ASTex::ImageRGB32;

template<typename I>
class Pcts
{
public:

    Pcts();

    void setTexture(const I& texture) {m_texture = texture;}

    void setWidth(int width) {m_width=width;}
    void setHeight(int height) {m_height=height;}
    void setSize(int width, int height) {m_width=width; m_height=height;}

    void setNbPasses(unsigned nbPasses) {m_nbPasses=nbPasses;}

    void setMinimumSizeLog(unsigned logValue) {m_minimumSizeLog=logValue;}

    void setLvl0BlockSize(unsigned blockSize) {m_lvl0BlockSize=blockSize;}

    void generate();

private:

    I m_texture;
    int m_width;
    int m_height;
    unsigned m_nbPasses;
    unsigned m_minimumSizeLog;
    unsigned m_lvl0BlockSize;
};

template<typename I>
Pcts<I>::Pcts() :
    m_texture(),
    m_width(800),
    m_height(800),
    m_nbPasses(2),
    m_minimumSizeLog(4)
{}

template<typename I>
void Pcts<I>::generate()
{
    Mipmap<I> pyramidInput;
    ImageIndex2 indexImageLevel0;

    //pyramid building

    unsigned maxReductionLevel = std::log2(std::min(m_texture.width(), m_texture.height())) - m_minimumSizeLog;

    pyramidInput.setTexture(m_texture);
    pyramidInput.setMode(ISOTROPIC);
    pyramidInput.setMaxPowReductionLevel(maxReductionLevel);
    pyramidInput.generate();

    const I& maxReductionMipmap = pyramidInput.mipmap(maxReductionLevel, maxReductionLevel);
    for(unsigned i=0; i<=maxReductionLevel; ++i)
        IO::save01_in_u8(pyramidInput.mipmap(i, i), std::string("/home/nlutz/pyramid") + std::to_string(i) + ".png");

    int blockSize=10;

    int     widthBlockyImage = m_width/std::pow(2.0, maxReductionLevel) + 1,
            heightBlockyImage= m_height/std::pow(2.0, maxReductionLevel) + 1;
    indexImageLevel0.initItk(widthBlockyImage, heightBlockyImage);
    for(int x=0; x<indexImageLevel0.width(); x+=blockSize)
        for(int y=0; y<indexImageLevel0.height(); y+=blockSize)
        {
            int xT=(rand())%(maxReductionMipmap.width()-blockSize-1);
            int yT=(rand())%(maxReductionMipmap.height()-blockSize-1);
            for(int x2=x; x2<blockSize+x && x2<indexImageLevel0.width(); ++x2)
                for(int y2=y; y2<blockSize+y && y2<indexImageLevel0.height(); ++y2) //TODO: preconditions
                {
                    indexImageLevel0.pixelAbsolute(x2, y2)[0] = xT+x2-x;
                    indexImageLevel0.pixelAbsolute(x2, y2)[1] = yT+y2-y;
                }
        }
    I imageLevel0;

    auto turnIndexIntoImage = [&] (const ImageIndex2& indexImage, I& image, int s)
    {
        image.initItk(indexImage.width(), indexImage.height());
        image.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
        {
            const I& mipmap = pyramidInput.mipmap(s, s);
            const ImageIndex2::PixelType &pixIndirect = indexImage.pixelAbsolute(x, y);
            pix = mipmap.pixelAbsolute(pixIndirect[0], pixIndirect[1]);
        });
    };

    turnIndexIntoImage(indexImageLevel0, imageLevel0, maxReductionLevel);
    IO::save01_in_u8(imageLevel0, "/home/nlutz/testBlock.png");

    //Synthesis
    for(int s=maxReductionLevel-1; s>=0; --s)
    {
        ImageIndex2 indexImageLevel1;
        indexImageLevel1.initItk(indexImageLevel0.width()*2, indexImageLevel0.height()*2);
        indexImageLevel0.for_all_pixels([&] (const typename I::PixelType &pix, int x, int y)
        {
            indexImageLevel1.pixelAbsolute(2*x, 2*y) = pix*2;

            indexImageLevel1.pixelAbsolute(2*x+1, 2*y) = pix*2;
            indexImageLevel1.pixelAbsolute(2*x+1, 2*y)[0]++;

            indexImageLevel1.pixelAbsolute(2*x, 2*y+1) = pix*2;
            indexImageLevel1.pixelAbsolute(2*x, 2*y+1)[1]++;

            indexImageLevel1.pixelAbsolute(2*x+1, 2*y+1) = pix*2;
            indexImageLevel1.pixelAbsolute(2*x+1, 2*y+1)[0]++;
            indexImageLevel1.pixelAbsolute(2*x+1, 2*y+1)[1]++;
        });
        I imageLevel1;
        imageLevel1.initItk(indexImageLevel1.width(), indexImageLevel1.height());
        turnIndexIntoImage(indexImageLevel1, imageLevel1, s);

        IO::save01_in_u8(imageLevel1, std::string("/home/nlutz/testBlock") + std::to_string(s) +".png");
        indexImageLevel0 = indexImageLevel1;
    }

}

}

#endif
