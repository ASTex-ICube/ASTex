#ifndef __PCTS_H__
#define __PCTS_H__

#include "ASTex/mipmap.h"
#include "ASTex/easy_io.h"
#include "ASTex/Stamping/sampler.h"

#define PCTS_DEBUG_DIRECTORY "/home/nlutz/"

namespace ASTex
{

using ImageIndex2 = ASTex::ImageRGB32;

double mse(const ImageRGBd& i1, const ImageRGBd& i2, int x1, int y1, int x2, int y2, int v)
{
    double error=0.0;
    static ImageRGBd::PixelType ms_zero;
    ImageRGBd::PixelType diff=ms_zero;

    unsigned hit=0;
    for(int dx=-v; dx<=v; ++dx)
        for(int dy=-v; dy<=v; ++dy)
        {
            int xx1 = x1+dx, xx2 = x2+dx;
            int yy1 = y1+dy, yy2 = y2+dy;

            if(xx1 >= 0 && xx1<i1.width() &&
               xx2 >= 0 && xx2<i2.width() &&
               yy1 >= 0 && yy1<i1.height() &&
               yy2 >= 0 && yy1<i2.height())
            {
                for(int i=0; i<3; ++i)
                {
                    diff[i] += ((i1.pixelAbsolute(xx1, yy1)[i]) - i2.pixelAbsolute(xx2, yy2)[i]) *
                                ((i1.pixelAbsolute(xx1, yy1)[i]) - i2.pixelAbsolute(xx2, yy2)[i]);
                }
                ++hit;
            }
        }
    assert(hit!=0 && "mse: error with parameters");
    return error=(diff[0] + diff[1] + diff[2])/3.0/hit;
}

template<typename I>
class Pcts
{
public:

    Pcts();
    Pcts(const I& texture);

    void setTexture(const I& texture) {m_texture = texture; m_textureSet=true;}

    void setWidth(int width) {m_width=width;}
    void setHeight(int height) {m_height=height;}
    void setSize(int width, int height) {m_width=width; m_height=height;}

    void setNbPasses(unsigned nbPasses) {m_nbPasses=nbPasses;}

    void setMinimumSizeLog(unsigned logValue) {m_minimumSizeLog=logValue;}

    void setLvl0BlockSize(unsigned blockSize) {m_lvl0BlockSize=blockSize;}

    void setCorrectionNeighborhood(unsigned neighborHood) {m_neighborhood=neighborHood;}

    void setNbSamplesNNM(unsigned nbSamples) {m_nbSamplesNNM=nbSamples;}
    void setRadiusScaleNNM(unsigned radiusScale) {m_radiusScaleNNM=radiusScale;}
    void setNbRefinementsNNM(unsigned nbRefinements) {m_nbRefinementsNNM=nbRefinements;}

    I generate();

private:

    I m_texture;
    int m_width;
    int m_height;
    unsigned m_nbPasses;
    unsigned m_minimumSizeLog;
    unsigned m_lvl0BlockSize;

    unsigned m_neighborhood;

    unsigned m_nbSamplesNNM;
    unsigned m_radiusScaleNNM;
    unsigned m_nbRefinementsNNM;

    bool m_textureSet;
};

template<typename I>
Pcts<I>::Pcts() :
    m_texture(),
    m_width(800),
    m_height(800),
    m_nbPasses(2),
    m_minimumSizeLog(4),
    m_lvl0BlockSize(10),
    m_neighborhood(2),
    m_nbSamplesNNM(5),
    m_radiusScaleNNM(8),
    m_nbRefinementsNNM(2),
    m_textureSet(false)
{}

template<typename I>
Pcts<I>::Pcts(const I& texture) :
    m_texture(texture),
    m_width(800),
    m_height(800),
    m_nbPasses(2),
    m_minimumSizeLog(4),
    m_lvl0BlockSize(10),
    m_neighborhood(2),
    m_nbSamplesNNM(5),
    m_radiusScaleNNM(8),
    m_nbRefinementsNNM(2),
    m_textureSet(true)
{}

template<typename I>
I Pcts<I>::generate()
{
    assert(m_textureSet &&
           "Pcts::generate: a texture must be set (try using Pcts::setTexture()).");

    Mipmap<I> pyramidInput;
    ImageIndex2 indexImageLevel0;

    //pyramid building

    unsigned maxReductionLevel = std::log2(std::min(m_texture.width(), m_texture.height())) - m_minimumSizeLog;

    pyramidInput.setTexture(m_texture);
    pyramidInput.setMode(ISOTROPIC);
    pyramidInput.setMaxPowReductionLevel(maxReductionLevel);
    pyramidInput.generate();

    const I& lvl0mipmap = pyramidInput.mipmap(maxReductionLevel, maxReductionLevel);
    for(unsigned i=0; i<=maxReductionLevel; ++i)
        IO::save01_in_u8(pyramidInput.mipmap(i, i), std::string(PCTS_DEBUG_DIRECTORY) + "pyramid" + std::to_string(i) + ".png");

    int     widthBlockyImage = m_width/std::pow(2.0, maxReductionLevel) + 1,
            heightBlockyImage= m_height/std::pow(2.0, maxReductionLevel) + 1;
    indexImageLevel0.initItk(widthBlockyImage, heightBlockyImage);
    for(int x=0; x<indexImageLevel0.width(); x+=m_lvl0BlockSize)
        for(int y=0; y<indexImageLevel0.height(); y+=m_lvl0BlockSize)
        {
            int xT=(rand())%(std::max(1, int(lvl0mipmap.width()-m_lvl0BlockSize)));
            int yT=(rand())%(std::max(1, int(lvl0mipmap.height()-m_lvl0BlockSize)));
            for(unsigned x2=x; x2<m_lvl0BlockSize+x && x2<(unsigned)indexImageLevel0.width(); ++x2)
                for(unsigned y2=y; y2<m_lvl0BlockSize+y && y2<(unsigned)indexImageLevel0.height(); ++y2) //TODO: preconditions
                {
                    indexImageLevel0.pixelAbsolute(x2, y2)[0] = xT+x2-x;
                    indexImageLevel0.pixelAbsolute(x2, y2)[1] = yT+y2-y;
                }
        }
    I imageLevel0;

    /**
      * Lambda lmbd_lookupIndexIntoImage performs the lookup pyramid.mipmap(s) o indexImage and copies the result to image.
    **/
    auto lmbd_lookupIndexIntoImage = [&] (const ImageIndex2& indexImage, I& image, int s)
    {
        image.initItk(indexImage.width(), indexImage.height());
        image.for_all_pixels([&] (typename I::PixelType &pix, int x, int y)
        {
            const I& mipmap = pyramidInput.mipmap(s, s);
            const ImageIndex2::PixelType &pixIndirect = indexImage.pixelAbsolute(x, y);
            pix = mipmap.pixelAbsolute(pixIndirect[0], pixIndirect[1]);
        });
        return;
    };

    /**
      * Lambda lmbd_fakeColorIndex turns image into a vizualizable index map,
      * of colors (0,0)=red, (1,0)=black, (0,1)=yellow, (1,1)=green.
    **/
    auto lmbd_fakeColorIndex = [&] (const ImageIndex2& indexImage, ImageRGBd& image)
    {
        image.initItk(indexImage.width(), indexImage.height());
        ImageRGBd::PixelType color0, color1, color2;
        color0[0]=1; color0[1]=0; color0[2]=0;
        color1[0]=0; color1[1]=1; color1[2]=0;
        color2[0]=1; color2[1]=1; color2[2]=0;
        float s, t;
        image.for_all_pixels([&] (ImageRGBd::PixelType &pix, int x, int y)
        {
            s=indexImage.pixelAbsolute(x, y)[0]/float(image.width());
            t=indexImage.pixelAbsolute(x, y)[1]/float(image.height());
            pix = color0*(1-s)*(1-t) + color1*s*t + color2*(1-s)*t;
        });
        return;
    };

    lmbd_lookupIndexIntoImage(indexImageLevel0, imageLevel0, maxReductionLevel);
    IO::save01_in_u8(imageLevel0, std::string(PCTS_DEBUG_DIRECTORY) + "testBlock" + std::to_string(maxReductionLevel) + ".png");

    //Synthesis
    for(int s=maxReductionLevel; s>=0; --s)
    {
        //Correction

        lmbd_lookupIndexIntoImage(indexImageLevel0, imageLevel0, s);

        for(unsigned nx=0; nx<m_neighborhood; ++nx)
            for(unsigned ny=0; ny<m_neighborhood; ++ny)
                for(unsigned x=nx; x<unsigned(indexImageLevel0.width()); x+=m_neighborhood)
                    for(unsigned y=ny; y<unsigned(indexImageLevel0.height()); y+=m_neighborhood)
                    {
                        itk::Index<2> idErrMin;
                        double errMin;

                        idErrMin[0]=indexImageLevel0.pixelAbsolute(x, y)[0];
                        idErrMin[1]=indexImageLevel0.pixelAbsolute(x, y)[1];

                        errMin=mse(pyramidInput.mipmap(s, s), imageLevel0, idErrMin[0], idErrMin[1], x, y, m_neighborhood);

                        /**
                          * Lambda updateErrorMin updates everything related to the MSE,
                          * including errMin, idErrMin, and the concerned textures.
                        **/
                        auto updateErrorMin = [&pyramidInput, &imageLevel0, &indexImageLevel0, &idErrMin, &errMin, s]
                                (double newErrMin, itk::Index<2> newIdErrMin, int x, int y)
                        {
                            errMin = newErrMin;
                            idErrMin = newIdErrMin;

                            indexImageLevel0.pixelAbsolute(x, y)[0]=idErrMin[0];
                            indexImageLevel0.pixelAbsolute(x, y)[1]=idErrMin[1];

                            imageLevel0.pixelAbsolute(x, y)=pyramidInput.mipmap(s, s).pixelAbsolute(idErrMin);
                        };

                        /**
                          * Lambda nearestNeighborMatch (NMM) searches for the nnm in x and y.
                          * dx and dy are used to easily expand the borders of each block.
                        **/
                        auto nearestNeighborMatch = [&] (int dx, int dy)
                        {
                            int px=indexImageLevel0.pixelAbsolute(x+dx, y+dy)[0];
                            int py=indexImageLevel0.pixelAbsolute(x+dx, y+dy)[1];

                            double err2=mse(pyramidInput.mipmap(s, s), imageLevel0, px, py, x, y, m_neighborhood);
                            if(err2<errMin)
                            {
                                itk::Index<2> pxpy;
                                pxpy[0]=px;
                                pxpy[1]=py;
                                updateErrorMin(err2, pxpy, x, y);
                            }
                            Stamping::SamplerPoisson sp;
                            sp.setNbPoints(m_nbSamplesNNM);
                            sp.setGenerateInCircle(true);
                            int radius;
                            for(unsigned n=0; n<m_nbRefinementsNNM &&
                                              (radius=m_radiusScaleNNM*m_neighborhood/pow(2.0, n))>0;
                                ++n)
                            {
                                std::vector<Eigen::Vector2f> spResult=sp.generate();

                                itk::Index<2> idPoisson;
                                for(unsigned i=0; i<spResult.size(); ++i)
                                {
                                    idPoisson[0]=(2*spResult[i][0] - 1)*radius + px;
                                    idPoisson[1]=(2*spResult[i][1] - 1)*radius + py;
                                    if(idPoisson[0]>=0 && idPoisson[0]<pyramidInput.mipmap(s, s).width() &&
                                       idPoisson[1]>=0 && idPoisson[1]<pyramidInput.mipmap(s, s).height())
                                    {
                                        err2=mse(pyramidInput.mipmap(s, s), imageLevel0, idPoisson[0], idPoisson[1],
                                                        x, y, m_neighborhood);
                                        if(err2<errMin)
                                            updateErrorMin(err2, idPoisson, x, y);
                                    }
                                }
                            }
                        };

                        nearestNeighborMatch(0, 0);
                        if(y+1<(unsigned)indexImageLevel0.height())
                            nearestNeighborMatch(0, 1);
                        if(x+1<(unsigned)indexImageLevel0.width())
                            nearestNeighborMatch(1, 0);
                        if(x+1<(unsigned)indexImageLevel0.width() && y+1<(unsigned)indexImageLevel0.height())
                        nearestNeighborMatch(1, 1);
                    }

        if(s>0)
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
            lmbd_lookupIndexIntoImage(indexImageLevel1, imageLevel1, s-1);
            IO::save01_in_u8(imageLevel1, std::string(PCTS_DEBUG_DIRECTORY) + "testBlock" + std::to_string(s) + ".png");

            ImageRGBd fakeColorIndexImageLevel0;
            lmbd_fakeColorIndex(indexImageLevel0, fakeColorIndexImageLevel0);

            indexImageLevel0 = indexImageLevel1;

            IO::save01_in_u8(fakeColorIndexImageLevel0, std::string(PCTS_DEBUG_DIRECTORY) + "fakeColorIndex" + std::to_string(s) + ".png");
        }

        ImageRGBd fakeColorIndexImageLevel0;
        lmbd_fakeColorIndex(indexImageLevel0, fakeColorIndexImageLevel0);

        IO::save01_in_u8(fakeColorIndexImageLevel0, std::string(PCTS_DEBUG_DIRECTORY) + "fakeColorIndex" + std::to_string(s) + ".png");
    }

    return imageLevel0;
}

}

#endif
