#ifndef __RESIDUAL_H__
#define __RESIDUAL_H__

#include "ASTex/image_gray.h"

namespace ASTex
{

namespace Sparse
{

class Residual
{
public:

    using ImageType = ImageGrayd;
    using SizeType = itk::Size;

    Residual();

    void setSize(SizeType size);
    void setNbLayers(unsigned nbLayers);
    ImageType &layer(unsigned layerId);
    const ImageType &layer(unsigned layerId) const;

private:

    void assert_integrity();

    std::vector<ImageGrayd> m_layers;
    SizeType m_size;
};

}

}

#endif
