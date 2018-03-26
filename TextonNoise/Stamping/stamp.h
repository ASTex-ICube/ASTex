
#ifndef __STAMP__H__
#define __STAMP__H__

#include "vector"

#include <vector>
#include <iostream>
#include <fstream>
#include <memory.h>

#include <cmath>
#include <ASTex/special_io.h>
#include <ASTex/easy_io.h>
#include <ASTex/fourier.h>
#include <ASTex/local_spectrum.h>
#include <ASTex/utils.h>
#include <ASTex/distances_maps.h>
#include <Eigen/Core>

namespace ASTex
{

namespace Stamping
{

//Declaration of StampBase
//TODO: for every class, constrain I to be an ImageCommon
template<typename I>
class StampBase
{
public:
    StampBase() {}
    StampBase(const StampBase &other) : StampBase() {}

    typedef typename I::PixelType PixelType;

    virtual PixelType   pixel(double x, double y) const = 0;
    virtual Size        size() const = 0;
    virtual int         width() const = 0;
    virtual int         height() const = 0;

private:
};

//Declaration of StampDiscrete
//TODO: add mask
template<typename I>
class StampDiscrete : public StampBase<I>
{
public:
    StampDiscrete();
    StampDiscrete(const StampDiscrete &other);
    StampDiscrete(const I &image);

    //types

    typedef typename StampBase<I>::PixelType PixelType;
    class Dimensions
    {
    public:
        Dimensions(double x=1.0, double y=1.0): dimX(x), dimY(y) {}
        double dimX;
        double dimY;
    };
    typedef enum{SD_BILINEAR, SD_NEAREST, SD_TRUNC} interpolation_rule_t;

    //get

    Size        size() const                        {return m_image.size();}
    int         width() const                       {return m_image.width();}
    int         height() const                      {return m_image.height();}

    Dimensions  dimensions() const                  {return m_dimensions;}
    double dimX() const                             {return m_dimensions.dimX;}
    double dimY() const                             {return m_dimensions.dimY;}
    interpolation_rule_t interpolationRule() const  {return m_rule;}

    const I& image() const                          {return m_image;}

    PixelType   pixel(double x, double y) const;

    //set

    void setInterpolationRule(interpolation_rule_t rule);
    void setDimensions(Dimensions dimensions);
    void setDimX(double dimX);
    void setDimY(double dimY);

    void setImage(const I& image); //TODO: might be sensitive, check for itk documentation

private:

    inline typename I::PixelType bilinear_interpolation_(typename I::PixelType q11,
                                                         typename I::PixelType q12,
                                                         typename I::PixelType q21,
                                                         typename I::PixelType q22,
                                                         int x1, int y1,
                                                         double x, double y);

    I                       m_image;
    interpolation_rule_t    m_rule;
    Dimensions              m_dimensions;

    static PixelType ms_zero;
};

//Declaration of StampContinuous?

//Implementation of StampDiscrete

template<typename I>
StampDiscrete<I>::StampDiscrete() :
    m_image(),
    m_rule(SD_BILINEAR),
    m_dimensions(Dimensions(1.0, 1.0))
{}

template<typename I>
StampDiscrete<I>::StampDiscrete(const StampDiscrete &other) :
    m_image(other.image()), //warning!!!
    m_rule(other.interpolationRule()),
    m_dimensions(other.dimensions())
{}

template<typename I>
typename StampDiscrete<I>::PixelType StampDiscrete<I>::ms_zero;

template<typename I>
StampDiscrete<I>::StampDiscrete(const I &image) :
    StampBase<I>::StampBase(),
    m_image(image), //TODO: copy pixels?
    m_rule(SD_BILINEAR),
    m_dimensions(Dimensions(1.0, 1.0))
{
//    m_image.initItk(image.width(), image.height());
//    m_image.copy_pixels(image);
}

template<typename I>
typename StampDiscrete<I>::PixelType StampDiscrete<I>::pixel(double x, double y) const
{
    int dx, dy, tx, ty;
    PixelType pixel = ms_zero;
    double dimX, dimY;
    typename I::PixelType q11, q12, q21, q22;

    dimX = x / m_dimensions.dimX;
    dimY = y / m_dimensions.dimY; //set x and y to scale;
    tx = (int) dimX;
    ty = (int) dimY;

    auto lmbd_is_in_range = [&] (int lx, int ly) -> bool
    {
        return lx >= 0 && lx < m_image.width() && ly >= 0 && ly < m_image.height();
    };

    if(m_rule == SD_BILINEAR)
    {
        q11 = tx >= 0 && ty >= 0 && tx < m_image.width() && ty < m_image.height() ? m_image.pixelAbsolute(tx, ty) : ms_zero;
        q12 = tx >= 0 && ty >= -1 && tx < m_image.width() && ty < m_image.height()-1 ? m_image.pixelAbsolute(tx, ty+1) : ms_zero;
        q21 = tx >= -1 && ty >= 0 && tx < m_image.width()-1 && ty < m_image.height() ? m_image.pixelAbsolute(tx+1, ty) : ms_zero;
        q22 = tx >= -1 && ty >= -1 && tx < m_image.width()-1 && ty < m_image.height()-1 ? m_image.pixelAbsolute(tx+1, ty+1) : ms_zero;

        pixel = bilinear_interpolation_(q11, q12, q21, q22, tx, ty, dimX, dimY);
    }
    else if(m_rule == SD_NEAREST) //TODO: that's not actually a real nearest is it
    {
        dx = dimX - tx;
        dy = dimY - ty;
        if(dx >= 0.5)
            ++tx;
        if(dy >= 0.5)
            ++ty;

        pixel = lmbd_is_in_range(tx, ty) ? m_image.pixelAbsolute(tx, ty) : ms_zero;
    }
    else //m_rule == SD_TRUNC
        pixel = lmbd_is_in_range(tx, ty) ? m_image.pixelAbsolute(tx, ty) : ms_zero;

    return pixel;
}

template<typename I>
typename I::PixelType StampDiscrete<I>::bilinear_interpolation_(typename I::PixelType q11,
                                                                typename I::PixelType q12,
                                                                typename I::PixelType q21,
                                                                typename I::PixelType q22,
                                                                int x1, int y1,
                                                                double x, double y)
{
    double x2x, y2y, yy1, xx1;
    yy1 = y - y1;
    xx1 = x - x1;
    x2x = 1.0 - xx1;
    y2y = 1.0 - yy1;
    return
        q11 * (x2x * y2y) +
        q21 * (xx1 * y2y) +
        q12 * (x2x * yy1) +
        q22 * (xx1 * yy1);
}

//set

template<typename I>
void StampDiscrete<I>::setInterpolationRule(StampDiscrete<I>::interpolation_rule_t rule)
{
    m_rule = rule;
}

template<typename I>
void StampDiscrete<I>::setDimensions(StampDiscrete<I>::Dimensions dimensions)
{
    assert(dimensions.dimX>0 && dimensions.dimY>0
           && "StampDiscrete::setDimensions: dimensions must be > 0");
    m_dimensions = dimensions;
}

template<typename I>
void StampDiscrete<I>::setDimX(double dimX)
{
    assert(dimX>0
           && "setDimX: dimX must be > 0");
    m_dimensions.dimX = dimX;
}

template<typename I>
void StampDiscrete<I>::setDimY(double dimY)
{
    assert(dimY>0
           && "setDimY: dimY must be > 0");
    m_dimensions.dimY = dimY;
}

template<typename I>
void StampDiscrete<I>::setImage(const I& image)
{
    assert(image.isInitialized());
    m_image = image;
}

} //namespace Stamping

} //namespace ASTex


#endif //__STAMP__H__
