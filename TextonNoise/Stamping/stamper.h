
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

class Stamper
{
    protected:
        std::vector<Eigen::Vector2f> m_pointArray;

        ImageRGBd m_stamp;

    public:

        Stamper(const std::vector<Eigen::Vector2f> &pointArray, const ImageRGBd &tampon);

        virtual ImageRGBd generate(int imageWidth, int imageHeight) = 0;
};

class BombingStamper : public Stamper
{
public:

    BombingStamper(const std::vector<Eigen::Vector2f> &pointArray, const ImageRGBd &tampon);

    ImageRGBd generate(int imageWidth, int imageHeight);
};

class TextonStamper : public Stamper
{
public:

    TextonStamper(const std::vector<Eigen::Vector2f> &pointArray, const ImageRGBd &tampon);

    ImageRGBd generate(int imageWidth, int imageHeight);

    //get

    double ratioX() {return m_ratioX;}
    double ratioY() {return m_ratioY;}
    bool periodicity() {return m_periodicity;}
    bool bilinearInterpolation() {return m_bilinearInterpolation;}
    bool useMargins() {return m_useMargins;}


    //set

    void setRatioX(double ratioX) {m_ratioX=ratioX; assert(ratioX>0 && "TextonStamper::setRatioX: ratioX must be > 0");}
    void setRatioY(double ratioY) {m_ratioY=ratioY; assert(ratioY>0 && "TextonStamper::setRatioY: ratioY must be > 0");}
    void setPeriodicity(bool periodicity) {m_periodicity = periodicity;}
    void setBilinearInterpolation(bool bi) {m_bilinearInterpolation = bi;}
    void setUseMargins(bool use) {m_useMargins = use;}

private:

    double m_ratioX; //< these define the resolution of the texton relative to the resolution of the texture (in X and Y)
    double m_ratioY;

    bool m_periodicity; //< if the stamping process is periodic or not.

    bool m_bilinearInterpolation; //< if the stamping process is allowed to stamp in between pixels or not.

    bool m_useMargins; //< if m_periodicity is true, defines if we don't extend the domain of the output with margins to improve the energy at the edges.
};

} //namespace Stamping

} //namespace ASTex


#endif //__STAMPER__H__
