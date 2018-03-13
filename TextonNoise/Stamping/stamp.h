
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

class StampBase
{
public:
    StampBase();

    virtual const ImageRGBd& sampledTexture()=0;
};

class StampDiscrete : public StampBase
{
public:
    StampDiscrete(const ImageRGBd &stamp);

    const ImageRGBd& sampledTexture();

private:

    ImageRGBd m_sampledTexture;
    ImageRGBd m_stamp;
};

} //namespace Stamping

} //namespace ASTex


#endif //__STAMP__H__
