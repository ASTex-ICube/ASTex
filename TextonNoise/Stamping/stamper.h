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

//TODO: consider removing the const of std::vector<const StampBase<I> *>
//if you think it should be authorized to modify the stamps given by the caller.

namespace ASTex
{

namespace Stamping
{

template<typename I>
class StamperBase
{
public:

    StamperBase(SamplerBase *sampler=0, const StampBase<I> *stamp=0);
    virtual ~StamperBase()                                      {delete(m_stampIndexGenerator);}

    //nested types

    typedef typename StampBase<I>::PixelType PixelType;
    class Func_stampIndexGeneratorBase
    {
    public:
        Func_stampIndexGeneratorBase() {}
        virtual ~Func_stampIndexGeneratorBase() {}
        virtual size_t operator() (void)=0;
    };
    class Func_stampIndexGeneratorZero : public Func_stampIndexGeneratorBase
    {
    public:
        Func_stampIndexGeneratorZero() : Func_stampIndexGeneratorBase() {}
        size_t operator() (void) {return 0;}
    };

    //main function

    virtual I generate(int imageWidth, int imageHeight) const = 0;

    //get

    SamplerBase* sampler() const  {return m_sampler;}
    const StampBase<I>* stamp(size_t index) const               {return m_stampVec[index];}

    const StampBase<I>*& operator[] (size_t index)              {return m_stampVec[index];}
    const StampBase<I>* const & operator[] (size_t index) const {return m_stampVec[index];}

    size_t stampVecSize() const                                 {return m_stampVec.size();}

    //set

    void setSampler(const SamplerBase* sampler)                 {m_sampler = sampler;}
    void setStamp(const StampBase<I>* stamp, size_t index)      {m_stampVec[index] = stamp;}

    //add/remove

    void addStamp(const StampBase<I>* stamp)                    {m_stampVec.push_back(stamp);}
    void removeStamp(size_t index)                              {m_stampVec.erase(m_stampVec.begin() + index); }

protected:

    void assert_null_before_generate_() const;

    SamplerBase                         *m_sampler;
    std::vector<const StampBase<I> *>   m_stampVec;
    Func_stampIndexGeneratorBase        *m_stampIndexGenerator;
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
StamperBase<I>::StamperBase(SamplerBase *sampler, const StampBase<I> *stamp)
    : m_sampler(sampler), m_stampVec(), m_stampIndexGenerator(0)
{
    if(stamp!=0) //allows the user to insert a default stamp
        m_stampVec.push_back(stamp);
}

template<typename I>
void StamperBase<I>::assert_null_before_generate_() const
{
    assert(m_sampler != NULL && "StamperBase::generate(w, h): sampler uninitialized (use setSampler(s) to give one)");
    for(typename std::vector<const StampBase<I> *>::const_iterator it=m_stampVec.begin(); it!=m_stampVec.end(); ++it)
        assert( (*it) != NULL && "StamperBase::generate(w, h): stamp uninitialized found (use addStamp(s) to give one)");
    return;
}

template<typename I>
StamperBombing<I>::StamperBombing(SamplerBase *sampler, const StampBase<I> *stamp) :
    StamperBase<I>(sampler, stamp)
{
    /// v it doesn't have to look this complicated if you define your own functor class for index generation.
    /// note that it will automatically be deleted by the base class destructor.
    this->m_stampIndexGenerator = new typename StamperBase<I>::Func_stampIndexGeneratorZero;
}

template<typename I>
I StamperBombing<I>::generate(int imageWidth, int imageHeight) const
{
    I im_out;
    int i, j, stampWidth, stampHeight, tx, ty, rx, ry;
    std::vector<Eigen::Vector2f> verticesArray; //stores the result of the sampler
    const StampBase<I> *stamp=0;
    this->assert_null_before_generate_(); /// < you can use this (or not) to check for a null sampler or a null stamp
    verticesArray = this->m_sampler->generate();

    im_out.initItk(imageWidth, imageHeight, true);

    /// v example of iteration over the sampler's result
    for(std::vector<Eigen::Vector2f>::const_iterator it=verticesArray.begin(); it!=verticesArray.end(); ++it)
    {
        stamp = this->m_stampVec[(*this->m_stampIndexGenerator)()]; /// < example of how to utilize a generator

        stampWidth = stamp->width();
        stampHeight = stamp->height();

        i = im_out.width() * (*it)[0]; //i & j: single point coordinates in im_out
        j = im_out.height() * (*it)[1];

        rx = i-stampWidth/2; //region origin coordinates (int)
        ry = j-stampHeight/2;
        Region reg = gen_region(rx, ry, stampWidth, stampHeight); /// < example of how to utilize regions for stamping
        im_out.for_region_pixels(reg, [&] (typename StampBase<I>::PixelType& pix, int x, int y)
        {
            if(x>0 && x<im_out.width() && y>0 && y<im_out.height())
            {
                tx=x-rx; //stamp coordinate shifted from the region origin coordinates
                ty=y-ry;
                pix = stamp->pixel(tx, ty); /// < example of mixing function (here, erase previous pixel with new one)
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
{
    this->m_stampIndexGenerator = new typename StamperBase<I>::Func_stampIndexGeneratorZero;
}

template<typename I>
I StamperTexton<I>::generate(int imageWidth, int imageHeight) const
{
    I im_out;
    std::vector<Eigen::Vector2f> verticesArray;
    double stampWidth, stampHeight;
    double i, j;
    double otx, oty; //texton origin in texture space (top left)
    double tx, ty; //texton coordinates
    //double dx, dy; //difference between int part of tx and double part of tx. Used to compute of the exact amount of energy which hit the texture
    double nbHit, nbHitPerPixel; //keeps track of the number of times the texture was hit
    double lambda; //lambda parameter in poisson process
    //assuming given texton is a zero mean texture with normalized variance

    const StampBase<I> *stamp=0;

    this->assert_null_before_generate_();
    verticesArray = this->m_sampler->generate();

    im_out.initItk(imageWidth, imageHeight, true);

    nbHit=0;

    for(std::vector<Eigen::Vector2f>::const_iterator it=verticesArray.begin(); it!=verticesArray.end(); ++it)
    {
        stamp = this->m_stampVec[(*this->m_stampIndexGenerator)()];
        stampWidth = stamp->width();
        stampHeight = stamp->height();

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

        nbHit += stampWidth*stampHeight; //with periodicity, the entire energy of the texton hits the texture all the time.
        //Without, we pretend it did and normalize including the size of the margins.
        if(m_periodicity)
        {
            im_out.for_region_pixels(reg, [&] (typename I::PixelType& pix, int x, int y) //with periodicity
            {
                (void)pix; /// pix is unused: remove parameter if there can be a (int, int) lambda
                tx=x-otx; //texton coordinate in texton space
                ty=y-oty; //texton coordinate

                im_out.pixelAbsolute((x+im_out.width())%imageWidth, (y+im_out.height())%imageHeight) += stamp->pixel(tx, ty);
            });
        }
        else
        {

            //the region we stamp : it is one pixel longer (per dim.) when the texton can stamped between pixels <=> when bilinear interpolation is activated
            //in othger terms, the region is a 2D bounding box for the texton
            im_out.for_region_pixels(reg, [&] (typename I::PixelType& pix, int x, int y) //without periodicity
            {
                (void)pix; /// pix is unused: remove parameter if there can be a (int, int) lambda
                if(x>=0 && y>=0 && x<imageWidth && y<imageHeight)
                { //here we compute pixel values, it's additive given there is 0 outside of the tx range.
                    tx=x-otx; //texton coordinate in texton space
                    ty=y-oty; //texton coordinate

                    im_out.pixelAbsolute(x, y) += stamp->pixel(tx, ty);

                }
            });
        }
    }
    //warning: assumes we use only one stamp / stamps with same width and height
    nbHitPerPixel = m_periodicity || !m_useMargins ? nbHit/(imageWidth*imageHeight) : nbHit/((imageWidth+stampWidth)*(imageHeight+stampHeight));

    lambda = double(nbHitPerPixel)/(stampWidth * stampHeight);

    std::cout << "Texton noise: number of impacts per pixel: " << std::to_string(nbHitPerPixel) << std::endl;

    im_out.for_all_pixels([&] (typename I::PixelType &pix) {
            pix = pix * (1.0/sqrt(lambda)) /*+ mean*/;
    });

    return im_out;
}

} //namespace Stamping

} //namespace ASTex


#endif //__STAMPER__H__
