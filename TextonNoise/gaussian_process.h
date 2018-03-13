#ifndef __GAUSSIAN_PROCESS_H_
#define __GAUSSIAN_PROCESS_H_

#include <ASTex/special_io.h>
#include <ASTex/fourier.h>
#include <ASTex/local_spectrum.h>
#include <ASTex/utils.h>
#include <ASTex/colorspace_filters.h>
#include <ASTex/mask.h>
#include <ASTex/distances_maps.h>
#include <ASTex/pca.h>
#include <ASTex/easy_io.h>

namespace ASTex
{

namespace Gaussian

{

//It's an ongoing project, I don't know how to make 2d gaussian process controlled by the PSD without FT+random phase

//template <typename REAL, typename std::enable_if<std::is_floating_point<REAL>::value>::type* = nullptr >
//void randomPhaseCosinus(const CommonSpectral<ImageGrayBase<REAL>, false>& modulus, ImageCommon<ImageGrayBase<REAL>, false> &output)
//{

//}

}

}

#endif

