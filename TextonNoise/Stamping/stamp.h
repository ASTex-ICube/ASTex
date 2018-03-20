
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
    StampBase();

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
    StampDiscrete(const I &stamp);

    //types

    typedef typename StampBase<I>::PixelType PixelType;
    class Dimensions
    {
        Dimensions(): dimX(1.0), dimY(1.0) {}
        double dimX;
        double dimY;
    };
    typedef enum{SD_BILINEAR, SD_NEAREST, SD_TRUNC} interpolation_rule_t;

    //get

    Size        size() const                        {return m_stamp.size();}
    int         width() const                       {return m_stamp.width();}
    int         height() const                      {return m_stamp.height();}

    Dimensions  dimensions() const                  {return m_dimensions;}
    double dimX() const                             {return m_dimensions.dimX;}
    double dimY() const                             {return m_dimensions.dimY;}
    interpolation_rule_t interpolationRule() const  {return m_rule;}

    PixelType   pixel(double x, double y) const;

    //set

    void setInterpolationRule(interpolation_rule_t rule);
    void setDimensions(Dimensions dimensions);
    void setDimX(double dimX);
    void setDimY(double dimY);

private:

    I                       m_stamp;
    interpolation_rule_t    m_rule;
    Dimensions              m_dimensions;
};

//Declaration of StampContinuous?

//Implementation of StampDiscrete

template<typename I>
StampDiscrete<I>::StampDiscrete() :
    m_stamp(),
    m_rule(SD_BILINEAR),
    m_dimensions(Dimensions(1.0, 1.0))
{}

template<typename I>
StampDiscrete<I>::StampDiscrete(const I &stamp) :
    m_stamp(stamp), //TODO: copy pixels?
    m_rule(SD_BILINEAR),
    m_dimensions(Dimensions(1.0, 1.0))
{}

template<typename I>
typename StampDiscrete<I>::PixelType StampDiscrete<I>::pixel(double x, double y) const
{
    int dx, dy, tx, ty;
    PixelType pixel;
    double dimX, dimY;

    dimX = x / m_dimensions.dimX;
    dimY = y / m_dimensions.dimY; //set x and y to scale;
    tx = (int) dimX;
    ty = (int) dimY;
    dx = dimX - tx;
    dy = dimY - ty;

    auto lmbd_is_in_range = [&] (int lx, int ly)
    {
        return lx >= 0 && lx < m_stamp.width() && ly >= 0 && ly < m_stamp.height();
    };

    if(m_rule == SD_BILINEAR)
    {
        if(tx >= 0 && ty >= 0)
            pixel += m_stamp.pixelAbsolute(tx, ty) * ((1-dx)*(1-dy));
        if(tx >=0 && ty < m_stamp.height()-1)
            pixel += m_stamp.pixelAbsolute(tx, ty+1) * ((1-dx)*dy);
        if(tx < m_stamp.width()-1 && ty >= 0)
            pixel += m_stamp.pixelAbsolute(tx+1, ty) * (dx*(1-dy));
        if(tx < m_stamp.width()-1 && ty < m_stamp.height()-1)
            pixel += m_stamp.pixelAbsolute(tx+1, ty+1) * (dx*dy);
    }
    else if(m_rule == SD_NEAREST)
    {
        if(dx >= 0.5)
            ++tx;
        if(dy >= 0.5)
            ++ty;

        pixel = lmbd_is_in_range(tx, ty) ? m_stamp.pixelAbsolute(tx, ty) : PixelType();
    }
    else //m_rule == SD_TRUNC
        pixel = lmbd_is_in_range(tx, ty) ? m_stamp.pixelAbsolute(tx, ty) : PixelType();

    return pixel;
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

} //namespace Stamping

} //namespace ASTex


#endif //__STAMP__H__
