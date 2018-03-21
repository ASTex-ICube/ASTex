
#ifndef __STAMPER__H__
#define __STAMPER__H__

#include <vector>

#include <vector>
#include <iostream>
#include <fstream>

#include <cmath>
#include <ASTex/easy_io.h>
#include <Eigen/Core>

#include <Stamping/stamp.h>
#include <Stamping/sampler.h>

namespace ASTex
{

namespace Stamping
{

template<typename I>
class StamperBase
{
public:

    StamperBase(SamplerBase *sampler=0, const StampBase<I> *stamp=0)
        : m_sampler(sampler), m_stamp(stamp) {}

    typedef typename StampBase<I>::PixelType PixelType;

    virtual I generate(int imageWidth, int imageHeight) const = 0;

    //get

    SamplerBase* sampler() const  {return m_sampler;}
    const StampBase<I>* stamp() const   {return m_stamp;}

    //set

    void setSampler(const SamplerBase* sampler)     {m_sampler = sampler;}
    void setStamp(const StampBase<I>* stamp)           {m_stamp = stamp;}

protected:

    void assert_null_before_generate_() const;

    SamplerBase                *m_sampler;
    const StampBase<I>         *m_stamp;

};

template<typename I>
class StamperBombing : public StamperBase<I>
{
public:

    StamperBombing(SamplerBase *sampler=0, const StampBase<I> *stamp=0);

    I generate(int imageWidth, int imageHeight) const;
};

template<typename I>
class StamperTexton : public StamperBase<I>
{
public:

    //warning: input must be a shifted and normalized texton (mean 0).
    StamperTexton(SamplerBase *sampler=0, const StampBase<I> *stamp=0);

    //warning: output is a shifted and normalized texton as well.
    I generate(int imageWidth, int imageHeight) const;

    //get

    double ratioX() {return m_ratioX;}
    double ratioY() {return m_ratioY;}
    bool periodicity() {return m_periodicity;}
    bool useMargins() {return m_useMargins;}


    //set

    void setRatioX(double ratioX) {m_ratioX=ratioX; assert(ratioX>0 && "TextonStamper::setRatioX: ratioX must be > 0");}
    void setRatioY(double ratioY) {m_ratioY=ratioY; assert(ratioY>0 && "TextonStamper::setRatioY: ratioY must be > 0");}
    void setPeriodicity(bool periodicity) {m_periodicity = periodicity;}
    void setUseMargins(bool use) {m_useMargins = use;}

private:

    double m_ratioX; //< these define the resolution of the texton relative to the resolution of the texture (in X and Y)
    double m_ratioY;

    bool m_periodicity; //< if the stamping process is periodic or not.

    bool m_useMargins; //< if m_periodicity is true, defines if we don't extend the domain of the output with margins to improve the energy at the edges.
};

template<typename I>
void StamperBase<I>::assert_null_before_generate_() const
{
    assert(m_sampler != NULL && "StamperBase::generate(w, h): sampler uninitialized (use setSampler(s) to give one)");
    assert(m_stamp != NULL && "StamperBase::generate(w, h): stamp uninitialized (use setStamp(s) to give one)");
    return;
}

template<typename I>
StamperBombing<I>::StamperBombing(SamplerBase *sampler, const StampBase<I> *stamp) :
    StamperBase<I>(sampler, stamp)
{}

template<typename I>
I StamperBombing<I>::generate(int imageWidth, int imageHeight) const
{
    I im_out;
    int i, j, textonWidth, textonHeight, tx, ty, rx, ry;
    this->assert_null_before_generate_();
    StampBase<I>::m_pointArray = StampBase<I>::m_sampler->generate();

    im_out.initItk(imageWidth, imageHeight, true);
    textonWidth = StampBase<I>::m_stamp->width();
    textonHeight = StampBase<I>::m_stamp->height();

    for(std::vector<Eigen::Vector2f>::const_iterator it=StampBase<I>::m_pointArray.begin(); it!=StampBase<I>::m_pointArray.end(); ++it)
    {
        i = im_out.width() * (*it)[0]; //i & j: single point coordinates in im_out
        j = im_out.height() * (*it)[1];

        rx = i-textonWidth/2; //region origin coordinates
        ry = j-textonHeight/2;
        Region reg = gen_region(rx, ry, textonWidth, textonHeight); //the region we stamp
        im_out.for_region_pixels(reg, [&] (typename StampBase<I>::PixelType& pix, int x, int y)
        {
            if(x>0 && x<im_out.width() && y>0 && y<im_out.height())
            {
                tx=x-rx; //texton coordinate shifted from the region origin coordinates
                ty=y-ry;
                pix = StampBase<I>::m_stamp->pixel(tx, ty);
            }
        });
    }

    return im_out;
}

template<typename I>
StamperTexton<I>::StamperTexton(SamplerBase *sampler, const StampBase<I> *stamp) :
    StamperBase<I>(sampler, stamp),
    m_ratioX(1.0),
    m_ratioY(1.0),
    m_periodicity(false),
    m_useMargins(true)
{}

template<typename I>
I StamperTexton<I>::generate(int imageWidth, int imageHeight) const
{
    ImageRGBd im_out;
    std::vector<Eigen::Vector2f> verticesArray;
    double stampWidth, stampHeight;
    double i, j;
    double otx, oty; //texton origin in texture space (top left)
    double tx, ty; //texton coordinates
    //double dx, dy; //difference between int part of tx and double part of tx. Used to compute of the exact amount of energy which hit the texture
    double nbHit, nbHitPerPixel; //keeps track of the number of times the texture was hit
    double lambda; //lambda parameter in poisson process
    //assuming given texton is a zero mean texture with normalized variance

    this->assert_null_before_generate_();
    verticesArray = this->m_sampler->generate();

    im_out.initItk(imageWidth, imageHeight, true);

    stampWidth = this->m_stamp->width();
    stampHeight = this->m_stamp->height();

    nbHit=0;

    for(std::vector<Eigen::Vector2f>::const_iterator it=verticesArray.begin(); it!=verticesArray.end(); ++it)
    {
        if(m_periodicity || (!m_periodicity && !m_useMargins))
        {
            i = imageWidth * (*it)[0];
            j = imageHeight * (*it)[1];
        }
        else
        {
            //we need to introduce margins here, to shoot textons outside of the domain
            i = (imageWidth + stampWidth ) * (*it)[0] - stampWidth/2.0;
            j = (imageHeight + stampHeight ) * (*it)[1] - stampHeight/2.0;
        }

        otx=i-stampWidth/2.0; //texton origin in texture space (top left)
        oty=j-stampHeight/2.0;

        Region reg = gen_region(std::floor(otx), std::floor(oty), stampWidth+1, stampHeight+1); //note: regions are weak when shooting between pixels

        nbHit += stampWidth*stampHeight; //with periodicity, the entire energy of the texton hits the texture all the time. Without, we pretend it did and normalize including the size of the margins.
        if(m_periodicity)
        {
            im_out.for_region_pixels(reg, [&] (ImageRGBd::PixelType& pix, int x, int y) //with periodicity
            {
                tx=x-otx; //texton coordinate in texton space
                ty=y-oty; //texton coordinate

                im_out.pixelAbsolute((x+im_out.width())%imageWidth, (y+im_out.height())%imageHeight) += this->m_stamp->pixel(tx, ty);
            });
        }
        else
        {

            //the region we stamp : it is one pixel longer (per dim.) when the texton can stamped between pixels <=> when bilinear interpolation is activated
            //in othger terms, the region is a 2D bounding box for the texton
            im_out.for_region_pixels(reg, [&] (ImageRGBd::PixelType& pix, int x, int y) //without periodicity
            {
                if(x>=0 && y>=0 && x<imageWidth && y<imageHeight)
                { //here we compute pixel values, it's additive given there is 0 outside of the tx range.
                    tx=x-otx; //texton coordinate in texton space
                    ty=y-oty; //texton coordinate

                    im_out.pixelAbsolute(x, y) += this->m_stamp->pixel(tx, ty);

                }
                /*

                if(x>=-stampWidth/2.0 && y>=-stampHeight/2.0 && x<imageWidth+stampWidth/2.0 && y<imageHeight+stampHeight/2.0)
                {//here we compute the energy, which is a good approximation of the mean number of impact. In the future, if a sampler gives the mean number of impact, use it!
                    dx=tx-std::floor(tx);
                    dy=ty-std::floor(ty);
                    if(tx<0) //then the rightmost part counts
                        if(ty<0) //then only the lower right part counts
                            nbHit += (1-dx)*(1-dy);
                        else if(ty>=stampHeight) //then only the upper right part counts
                            nbHit += (1-dx)*dy);
                        else //then only the rightmost part counts
                            nbHit += 1-dx;
                    else if(tx>=stampWidth) //then only the leftmost part counts
                        if(ty<0) //then only the lower left part counts
                            nbHit += dx*(1-dy);
                        else if(ty>=stampWidth) //then only the upper left part counts
                            nbHit += dx*dy;
                        else //then both the upper left and lower left parts count
                            nbHit += dx;
                    else if(ty<0) //then only the lower part
                        nbHit +=1-dy;
                    else if(ty>=stampHeight)


                } */
            });

        }
    }
    nbHitPerPixel = m_periodicity ? nbHit/(imageWidth*imageHeight) : nbHit/((imageWidth+stampWidth)*(imageHeight+stampHeight));

    lambda = double(nbHitPerPixel)/(this->m_stamp->width() * this->m_stamp->height());

    std::cout << "Texton noise: number of impacts per pixel: " << std::to_string(nbHitPerPixel) << std::endl;

    im_out.for_all_pixels([&] (ImageRGBd::PixelType &pix) {
        for(int i=0; i<3; ++i)
            pix[i] = 1.0/sqrt(lambda) * (pix[i]) /*+ mean[i]*/;
    });

    return im_out;
}

} //namespace Stamping

} //namespace ASTex


#endif //__STAMPER__H__
